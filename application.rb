require "./pubchem.rb"
require 'rack/session/dalli'
module OpenTox
  class Application < Service
    set :static, true
    set :root, File.dirname(__FILE__)
    also_reload './pubchem.rb'
    #enable :sessions
    use Rack::Session::Dalli, :cache => Dalli::Client.new

    @@pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"

    helpers do

=begin
      def pubchem_search url
        attempts = 0
        result = nil
        begin
          attempts += 1
          json = RestClient.get url, :timeout => 90000000
          result = JSON.parse json
          while result["Waiting"] do
            sleep 2
            listkey = result["Waiting"]["ListKey"]
            result = pubchem_search File.join(@pug_uri, "compound", "listkey", listkey, "cids", "JSON")
          end
        rescue
          if $!.message =~ /Timeout/i and attempts < 4
            sleep 2
            retry
          elsif $!.message =~ /Timeout/i and attempts >= 4
            File.open("timeouts","a+"){|f| f.puts url}
            puts url
            puts $!.message
          elsif $!.message.match /404/
            #not_found_error #TODO
          else
            puts url
            puts $!.message
          end
        end
      end
=end

      def image_uri cid
        File.join @@pug_uri, "compound", "cid", cid, "PNG"#?record_type=3d&image_size=small"
      end

    end

    before '/cid/:cid/*' do
      session[:compound] = PubChemCompound.new params[:cid] unless session[:compound] and session[:compound].cid == params[:cid]
    end

    get '/?' do
      haml :index
    end

    get '/cid/:cid/?' do
      session[:compound] = PubChemCompound.new params[:cid]
      haml :compound
    end

    get '/search/?' do
      #begin
        cids = RestClientWrapper.get(File.join(@@pug_uri,"compound","name",URI.escape(params[:name]),"cids","TXT")).split("\n")
        if cids.size == 1
          session[:compound] = PubChemCompound.new cids.first
          haml :compound
        elsif cids.size > 1
          @compounds = cids.collect{|cid| PubChemCompound.new cid }
          haml :select
        end
      #rescue
        #haml :not_found
      #end
    end

    get '/cid/:cid/targets/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].targets
      else
        @assays = PubChemCompound.new(params[:cid]).targets
      end
      haml :targets, :layout => false
    end

    get '/cid/:cid/nontargets/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].non_targets
      else
        @assays = PubChemCompound.new(params[:cid]).non_targets
      end
      haml :targets, :layout => false
    end

    get '/cid/:cid/other_active_assays/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].active_assays - session[:compound].targets
      else
        compound = PubChemCompound.new(params[:cid])
        @assays = compound.active_assays - compound.targets
      end
      haml :assays, :layout => false
    end

    get '/cid/:cid/other_inactive_assays/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].inactive_assays - session[:compound].non_targets
      else
        compound = PubChemCompound.new(params[:cid])
        @assays = compound.inactive_assays - compound.non_targets
      end
      haml :assays, :layout => false
    end

    get '/cid/:cid/predicted_targets/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].predicted_targets
      else
        @assays = PubChemCompound.new(params[:cid]).predicted_targets
      end
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/predicted_nontargets/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].predicted_non_targets
      else
        @assays = PubChemCompound.new(params[:cid]).predicted_non_targets
      end
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/other_predicted_active_assays/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].predicted_active_assays - session[:compound].predicted_targets
      else
        compound = PubChemCompound.new(params[:cid])
        @assays = compound.predicted_active_assays - compound.predicted_targets
      end
      haml :predicted_assays, :layout => false
    end

    get '/cid/:cid/other_predicted_inactive_assays/?' do
      if params[:cid] == session[:compound].cid
        @assays = session[:compound].predicted_inactive_assays - session[:compound].predicted_non_targets
      else
        compound = PubChemCompound.new(params[:cid])
        @assays = compound.predicted_inactive_assays - compound.predicted_non_targets
      end
      haml :assays, :layout => false
    end

    get '/cid/:cid/neighbors/?' do
      haml :neighbors, :layout => false
    end


    get '/cid/:cid/cosine/:cid2/?' do
      session[:compound].cosine(PubChemCompound.new(params[:cid2])).to_s
    end

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
