require "./pug.rb"
require "./pubchem.rb"
module OpenTox
  class Application < Sinatra::Base

    set :static, true
    set :root, File.dirname(__FILE__)

    configure :development do
      register Sinatra::Reloader
      also_reload './pubchem.rb'
    end

    before '/cid/:cid/*' do
      @compound = PubChemCompound.new params[:cid] 
    end

    get '/?' do
      haml :index
    end

    get '/cid/:cid/?' do
      haml :compound
    end

    get '/search/?' do
      @compounds = PubChemCompound.from_name params[:name]
      if @compounds.nil?
        haml :not_found
      elsif @compounds.is_a? Array
        haml :select
      else
        @compound = @compounds
        haml :compound
      end
    end

    get '/cid/:cid/targets/?' do
      @assays = @compound.targets
      if @assays.empty?
        "<br><em>No PubChem data</em></br>"
      else
        haml :targets, :layout => false
      end
    end

    get '/cid/:cid/nontargets/?' do
      @assays = @compound.non_targets
      if @assays.empty?
        "<br><em>No PubChem data</em></br>"
      else
        haml :targets, :layout => false
      end
    end

    get '/cid/:cid/other_active_assays/?' do
      @assays = @compound.active_assays - @compound.targets
      if @assays.empty?
        "<br><em>No PubChem data</em></br>"
      else
        haml :assays, :layout => false
      end
    end

    get '/cid/:cid/other_inactive_assays/?' do
      @assays = @compound.inactive_assays - @compound.non_targets
      if @assays.empty?
        "<br><em>No PubChem data</em></br>"
      else
        haml :assays, :layout => false
      end
    end

    get '/cid/:cid/predicted_targets/?' do
      @assays = @compound.predicted_targets
      puts @assays.inspect
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/predicted_nontargets/?' do
      @assays = @compound.predicted_non_targets
      haml :predicted_targets, :layout => false
    end

    get '/cid/:cid/other_predicted_active_assays/?' do
      @assays = @compound.predicted_active_assays - @compound.predicted_targets
      haml :predicted_assays, :layout => false
    end

    get '/cid/:cid/other_predicted_inactive_assays/?' do
      @assays = @compound.predicted_inactive_assays - @compound.predicted_non_targets
      haml :predicted_assays, :layout => false
    end

    get '/cid/:cid/neighbors/?' do
      haml :neighbors, :layout => false
    end

    get '/cid/:cid/cosine/:cid2/?' do
      @compound.cosine(PubChemCompound.new(params[:cid2])).round(3).to_s
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
