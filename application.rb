require 'base64'
require 'sinatra/base'
require "sinatra/reloader" 
require "rest-client"
require 'memcache'
require "yajl"
require 'yajl/json_gem'

class Application < Sinatra::Base

  set :static, true
  set :root, File.dirname(__FILE__)

  # doc @ http://pubchem.ncbi.nlm.nih.gov/pug_rest/
  PUG_URI = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
  SIMILARITY_THRESHOLD = 90
  MAX_NEIGHBORS = 100
  CACHE = MemCache.new 'localhost:11211'

  configure :development do
    register Sinatra::Reloader
  end

  helpers do

    def local route
      status, headers, body = call env.merge("PATH_INFO" => route)
      Yajl::Parser.parse body[0]
    rescue
      body[0] if body
    end

    def pubchem_search url
      attempts = 0
      begin
        attempts += 1
        puts url
        json = RestClient.get url, :timeout => 90000000
        Yajl::Parser.parse json
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

    def from_name name
      cids = local("/pug/name/#{CGI.escape(name)}")
    end

    [:fingerprint,:neighbors,:experiments,:predictions,:name].each do |method|
      send :define_method, method do |cid|
        local("/pug/cid/#{cid}/#{method}")
      end
    end

    def assays cid, outcome
      experiments(cid).select{|a| a["Activity Outcome"] == outcome} 
    rescue
    end

    def targets cid, outcome
      assays(cid, outcome).select{|a| a["Target GI"]}
    rescue
    end

    def predicted_assays cid, outcome
      case outcome
      when "active"
        predictions(cid).select{|a| a["p_active"] > a["p_inactive"]} 
      when "inactive"
        predictions(cid).select{|a| a["p_active"] < a["p_inactive"]} 
      end
    rescue
    end

    def predicted_targets cid, outcome
      predicted_assays(cid, outcome).select{|a| a["Target GI"]}
    rescue
    end

    def similarity cid1, cid2
      local("/pug/cid/#{cid1}/cosine/#{cid2}")
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

    def image_uri cid
      "/pug/cid/#{cid}/image"
    end
  end

  before '/pug/*' do
    content_type 'application/json'
    @result = CACHE.get request.path
    halt 200, @result unless @result.nil? # should be 304, but this does not work with local()
  end

  after '/pug/*' do
    CACHE.add request.path, @result, 72000 if @result
  end

  before '/cid/:cid/*' do
    @cid = params[:cid] 
  end

  get '/?' do
    haml :index
  end

  get '/cid/:cid/?' do
    @cid = params[:cid] 
    #local("/pug/cid/#{@cid}/predictions")
    pid = Process.fork{ local("/pug/cid/#{@cid}/predictions") }
    Process.detach pid
    haml :compound
  end

  get '/search/?' do
    redirect to("/") if params[:name].empty?
    @cids = from_name params[:name]
    if !@cids or @cids.empty?
      haml :not_found
    elsif @cids.size == 1
      redirect to("/cid/#{@cids.first}")
      #@cid = @cids.first
      #haml :compound
    else
      haml :select
    end
  end

  get '/cid/:cid/targets/:outcome' do
    @assays = targets params[:cid], params[:outcome]
    @assays and !@assays.empty? ? haml(:targets, :layout => false) : "<p><em>No PubChem data</em></p>" 
  end

  get '/cid/:cid/assays/:outcome' do
    @assays = assays(params[:cid], params[:outcome]) - targets(params[:cid], params[:outcome])
    @assays and !@assays.empty? ? haml(:assays, :layout => false) : "<p><em>No PubChem data</em></p>"
  end

  get '/cid/:cid/prediction/assays/:outcome' do
    @assays = predicted_assays(params[:cid], params[:outcome]) - predicted_targets(params[:cid], params[:outcome])
    @assays and !@assays.empty? ? haml(:predicted_assays, :layout => false) : "<p><em>Insuffucient PubChem data for read across predictions.</em></p>"
  end

  get '/cid/:cid/prediction/targets/:outcome' do
    @assays = predicted_targets params[:cid], params[:outcome]
    @assays and !@assays.empty? ? haml(:predicted_targets, :layout => false) : "<p><em>Insuffucient PubChem data for read across predictions.</em></p>"
  end

  get '/cid/:cid/neighbors/?' do
    haml :neighbors, :layout => false
  end

  get '/pug/cid/:cid/name' do
    @result = RestClient.get(File.join(PUG_URI, "compound", "cid", params[:cid], "property", "IUPACName","TXT")).chomp.to_json
  end

  get '/pug/cid/:cid/fingerprint' do
    # ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
    # it seems that only SDF formats contain fingerprints 
    sdf_lines = RestClient.get(File.join(PUG_URI, "compound", "cid", params[:cid], "SDF")).split("\n")
    index = sdf_lines.index(sdf_lines.grep(/PUBCHEM_CACTVS_SUBSKEYS/).first)
    @result = Base64.decode64(sdf_lines[index+1])[4..-1].unpack("B*").first[0..-8].split(//).collect{|c| c == "1"}.to_json
  end

  get '/pug/cid/:cid/experiments' do
    assays = []
    result = pubchem_search File.join(PUG_URI, "compound", "cid", params[:cid], "assaysummary", "JSON")
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

  get '/pug/cid/:cid/neighbors' do
    result = pubchem_search File.join(PUG_URI, "compound", "similarity", "cid", params[:cid], "JSON")+"?Threshold=#{SIMILARITY_THRESHOLD}&MaxRecords=#{MAX_NEIGHBORS}"
    while result["Waiting"] do
      sleep 2
      listkey = result["Waiting"]["ListKey"]
      result = pubchem_search File.join(PUG_URI, "compound", "listkey", listkey, "cids", "JSON")
    end
    result["IdentifierList"]["CID"].delete params[:cid].to_i
    #result["IdentifierList"]["CID"].each do |cid|
    #  @result << cid unless assays(cid,"active").empty? and assays(cid,"inactive").empty?
    #end
    @result = result["IdentifierList"]["CID"].to_json
  end

  get '/pug/name/:name' do
  begin
    @result = RestClient.get(File.join(PUG_URI,"compound","name",CGI.escape(params[:name]),"cids","TXT")).split("\n").to_json
  rescue
    @result = nil
  end
  end

  get '/pug/cid/:cid/image' do
    @result = RestClient.get File.join(PUG_URI, "compound", "cid", params[:cid], "PNG")
  end

  get '/pug/cid/:cid1/cosine/:cid2' do
    fp1 = fingerprint params[:cid1]
    fp2 = fingerprint params[:cid2]
    @result = cosine(fp1, fp2).to_json
  end

  get '/pug/cid/:cid/predictions' do
    assays = {}
    assay_details = {}
    neighbors(params[:cid]).each do |cid|
      neighbor_assays = experiments cid
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
#    RestClient.get File.join(PUG_URI,path), params
#  end

  get '/fp/?' do
    @fp = []
    YAML.load_file("false_positives.yaml").each do |pred|
      pred[:fp_targets].each do |gi,t|
        @fp << {
          "CID" => pred[:cid],
          "Target GI" => gi,
          "p_active" => t[:p][:active].first, 
          "p_inactive" => t[:p][:inactive].first, 
          :assays => t[:measured],
          :neighbors => t[:neighbors]
        }
      end
    end
    @fp.sort!{|a,b| b["p_active"] <=> a["p_active"]}
    haml :fp
  end

  # SASS stylesheet
  get '/style.css' do
    headers 'Content-Type' => 'text/css; charset=utf-8'
    scss :style
  end

end
