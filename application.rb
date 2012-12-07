require "./pubchem.rb"
require 'rack/session/dalli'
require "rack/cache"
module OpenTox
  class Application < Service
    set :static, true
    set :root, File.dirname(__FILE__)
    also_reload './pubchem.rb'

    @@pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"

    before do
      #cache_control :public, :max_age => 3600
    end

    before '/cid/:cid/*' do
      #cache_control :public, :max_age => 3600
      session[:compound] = PubChemCompound.new params[:cid] unless session[:compound] and session[:compound].cid == params[:cid]
    end

    get '/?' do
      #cache_control :public, :no_cache
      haml :index
    end

    get '/cid/:cid/?' do
      haml :compound
    end

    get '/search/?' do
      #cache_control :public, :no_cache
      begin
        cids = RestClient.get(File.join(@@pug_uri,"compound","name",CGI.escape(params[:name]),"cids","TXT")).split("\n")
        if cids.size == 1
          session[:compound] = PubChemCompound.new cids.first
          haml :compound
        elsif cids.size > 1
          @compounds = cids.collect{|cid| PubChemCompound.new cid }
          haml :select
        end
      rescue
        haml :not_found
      end
    end

    get '/cid/:cid/targets/?' do
      @assays = session[:compound].targets
      haml :targets, :layout => false
    end

    get '/cid/:cid/nontargets/?' do
      @assays = session[:compound].non_targets
      haml :targets, :layout => false
    end

    get '/cid/:cid/other_active_assays/?' do
      @assays = session[:compound].active_assays - session[:compound].targets
      haml :assays, :layout => false
    end

    get '/cid/:cid/other_inactive_assays/?' do
      @assays = session[:compound].inactive_assays - session[:compound].non_targets
      haml :assays, :layout => false
    end

    get '/cid/:cid/predicted_targets/?' do
      @assays = session[:compound].predicted_targets
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/predicted_nontargets/?' do
      @assays = session[:compound].predicted_non_targets
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/other_predicted_active_assays/?' do
      @assays = session[:compound].predicted_active_assays - session[:compound].predicted_targets
      haml :predicted_assays, :layout => false
    end

    get '/cid/:cid/other_predicted_inactive_assays/?' do
      @assays = session[:compound].predicted_inactive_assays - session[:compound].predicted_non_targets
      haml :predicted_assays, :layout => false
    end

    get '/cid/:cid/neighbors/?' do
      haml :neighbors, :layout => false
    end

    get '/cid/:cid/cosine/:cid2/?' do
      session[:compound].cosine(PubChemCompound.new(params[:cid2])).round(3).to_s
    end

=begin
    get '/aid/:aid/?' do
      puts File.join(@@pug_uri, "assay", "aid", params[:aid].to_s, "description", "JSON")
      json = RestClient.get File.join(@@pug_uri, "assay", "aid", params[:aid].to_s, "description", "JSON")
      @description = JSON.parse(json)["PC_AssayContainer"][0]["assay"]["descr"]
      haml :assay_description, :layout => false
    end

    get '/pubchem_proxy/*' do |path|
      puts path.inspect
      puts "http://pubchem.ncbi.nlm.nih.gov/rest/pug/#{path}"
      RestClientWrapper.get "http://pubchem.ncbi.nlm.nih.gov/rest/pug/#{path}"
    end
=end

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
  end
end
