require "json"
require 'base64'
require 'sinatra/base'
require "sinatra/reloader" 
require "rest-client"
require 'memcache'

class Application < Sinatra::Base

  # doc @ http://pubchem.ncbi.nlm.nih.gov/pug_rest/
  @@pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
  @@similarity_threshold = 90
  @@max_neighbors = 100
  
  CACHE = MemCache.new 'localhost:11211'

  helpers do

    def local route
      status, headers, body = call env.merge("PATH_INFO" => route)
      begin
        Float body[0]
      rescue
        JSON.parse body[0]
      end
    end

    def pubchem_search url
      attempts = 0
      begin
        attempts += 1
        puts url
        json = RestClient.get url, :timeout => 90000000
        JSON.parse json
      rescue
        if $!.message =~ /Timeout/i and attempts < 4
          sleep 2
          retry
        elsif $!.message =~ /Timeout/i and attempts >= 4
          File.open("timeouts","a+"){|f| f.puts url}
          puts url
          puts $!.message
          nil
        elsif $!.message.match /404/
          nil
        else
          puts url
          puts $!.message
          nil
        end
      end
    end

    def fingerprint cid
      local("/cid/#{cid}/fingerprint")
    end

    def neighbors cid
      local "/cid/#{cid}/neighbors"
    end

    def assays cid
      local "/cid/#{cid}/assays"
    end

    def similarity cid1, cid2
      local("/cid/#{cid1}/cosine/#{cid2}")
    end

    def cosine fp1, fp2
      if fp1 and fp2
        m11 = 0.0
        m01 = 0.0
        m10 = 0.0
        m00 = 0.0
        fp1.each_index do |i|
          m11 += 1 if (fp1[i] and fp2[i])
          m01 += 1 if (!fp1[i] and fp2[i]) 
          m10 += 1 if (fp1[i] and !fp2[i]) 
          m00 += 1 if (!fp1[i] and !fp2[i]) 
        end
        m11/((m01+m11)*(m10+m11))**0.5
      end
    end
  end

  before do
    @result = CACHE.get request.path
    halt 200, @result unless @result.nil? # should be 304, but this does not work with local()
  end

  after do
    CACHE.add request.path, @result, 7200
  end

  get '/cid/:cid/name' do
    @result = RestClient.get(File.join(@@pug_uri, "compound", "cid", params[:cid], "property", "IUPACName","TXT")).chomp
  end

  get '/cid/:cid/fingerprint' do
    # ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
    # it seems that only SDF formats contain fingerprints 
    sdf_lines = RestClient.get(File.join(@@pug_uri, "compound", "cid", params[:cid], "SDF")).split("\n")
    index = sdf_lines.index(sdf_lines.grep(/PUBCHEM_CACTVS_SUBSKEYS/).first)
    @result = Base64.decode64(sdf_lines[index+1])[4..-1].unpack("B*").first[0..-8].split(//).collect{|c| c == "1"}.to_json
  end

  get '/cid/:cid/assays' do
    assays = []
    result = pubchem_search File.join(@@pug_uri, "compound", "cid", params[:cid], "assaysummary", "JSON")
    if result and result["Table"]
      columns = result["Table"]["Columns"]["Column"]
      result["Table"]["Row"].collect{|cell| cell.values.flatten}.each do |row|
        assay = {}
        row.each_with_index do |cell,i|
          assay[columns[i]] = cell unless cell.empty? or columns[i] == "CID" 
        end
        assays << assay unless assay.empty?
      end
    end
    @result = assays.to_json
  end

  get '/cid/:cid/neighbors' do
    result = pubchem_search File.join(@@pug_uri, "compound", "similarity", "cid", params[:cid], "JSON")+"?Threshold=#{@@similarity_threshold}&MaxRecords=#{@@max_neighbors}"
    while result["Waiting"] do
      sleep 2
      listkey = result["Waiting"]["ListKey"]
      result = pubchem_search File.join(@@pug_uri, "compound", "listkey", listkey, "cids", "JSON")
    end
    result["IdentifierList"]["CID"].delete params[:cid].to_i
    @result = result["IdentifierList"]["CID"].to_json
  end

  get '/name/:name' do
    @result = RestClient.get(File.join(@@pug_uri,"compound","name",CGI.escape(params[:name]),"cids","TXT")).split("\n").to_json
  end

  get '/cid/:cid/image' do
    @result = RestClient.get File.join(@@pug_uri, "compound", "cid", params[:cid], "PNG")
  end

  get '/cid/:cid1/cosine/:cid2' do
    fp1 = fingerprint params[:cid1]
    fp2 = fingerprint params[:cid2]
    @result = cosine(fp1, fp2).to_json
  end

  get '/cid/:cid/predictions' do
    assays = {}
    assay_details = {}
    neighbors(params[:cid]).each do |cid|
      neighbor_assays = assays cid
      unless neighbor_assays.empty?
        neighbor_assays.each do |assay|
          if assay["Activity Outcome"] == "active" or assay["Activity Outcome"] == "inactive"
            assays[assay["AID"]] ||= []
            assays[assay["AID"]] << [cid,similarity(params[:cid],cid),assay["Activity Outcome"]]
            assay_details[assay["AID"]] ||= {}
            ["Target GI", "Target Name", "Assay Name"].each do |d|
              assay_details[assay["AID"]][d] = assay[d] #if assay[d]
            end
          end
        end
      end
    end
    predictions = []
    assays.each do |aid,neighbors|
      prediction = {"AID" => aid}
      neighbors.each do |neighbor|
        cid = neighbor[0]
        sim = neighbor[1]
        activity = neighbor[2]
        if activity == "active"
          prediction[:p_active] ? prediction[:p_active] = prediction[:p_active]*sim : prediction[:p_active] = sim
          prediction[:p_inactive] ? prediction[:p_inactive] = prediction[:p_inactive]*(1-sim) : prediction[:p_inactive] = 1-sim
        elsif activity == "inactive"
          prediction[:p_active] ? prediction[:p_active] = prediction[:p_active]*(1-sim) : prediction[:p_active] = 1-sim
          prediction[:p_inactive] ? prediction[:p_inactive] = prediction[:p_inactive]*sim : prediction[:p_inactive] = sim
        end
        ["Target GI", "Target Name", "Assay Name"].each do |d|
          prediction[d] = assay_details[aid][d] if assay_details[aid][d]
        end
      end
      predictions << prediction
    end
    @result = predictions.to_json
  end

#  get '/*' do |path|
#    RestClient.get File.join(@@pug_uri,path), params
#  end
end
