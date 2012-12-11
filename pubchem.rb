require '../opentox-client/lib/opentox-client.rb'
require 'json'
require 'base64'

def Math.gauss(x, sigma = 0.3) 
  d = 1.0 - x.to_f
  Math.exp(-(d*d)/(2*sigma*sigma))
end

module OpenTox

  class PubChemCompound < Compound
 
    attr_accessor :cid
    @@pug_proxy = "http://localhost:8081/"
    
    def initialize cid
      @cid = cid.to_s
    end

    def fingerprint
      JSON.parse RestClient.get(File.join(@@pug_proxy,"cid",@cid,"fingerprint"))
    end

    def self.from_name name
      cids = JSON.parse(RestClient.get(File.join(@@pug_proxy,"name",CGI.escape(name))))
      if cids.size == 1
        PubChemCompound.new cids.first
      elsif cids.empty?
        nil
      else
        cids.collect{|cid| PubChemCompound.new cid}
      end
    end

    def name
      RestClient.get(File.join(@@pug_proxy,"cid",cid,"name")).chomp.sub(/^"/,'').sub(/"$/,'')
    end

    def neighbors
      JSON.parse(RestClient.get(File.join(@@pug_proxy,"cid",@cid,"neighbors"))).collect{|n| PubChemCompound.new(n) }
    end

    def assays
      JSON.parse RestClient.get(File.join(@@pug_proxy,"cid",cid,"assays"))
    end

    def active_assays
      assays.select{|a| a["Activity Outcome"] == "active"} if assays
    end

    def inactive_assays
      assays.select{|a| a["Activity Outcome"] == "inactive"} if assays
    end

    def targets
      active_assays.select{|a| a["Target GI"]} if assays
    end

    def non_targets
      inactive_assays.select{|a| a["Target GI"]} if assays
    end

    def predicted_assays
      JSON.parse RestClient.get(File.join(@@pug_proxy,"cid",cid,"predictions"))
    end

    def predicted_active_assays
      predicted_assays.select{|a| a["p_active"] > a["p_inactive"]} if predicted_assays
    end

    def predicted_inactive_assays
      predicted_assays.select{|a| a["p_active"] < a["p_inactive"]} if predicted_assays
    end

    def predicted_targets
      predicted_active_assays.select{|a| a["Target GI"]} if predicted_assays
    end

    def predicted_non_targets
      predicted_inactive_assays.select{|a| a["Target GI"]} if predicted_assays
    end

    def image_uri
      File.join @@pug_proxy, "cid", @cid, "image"
    end

    def similarity compound
      cosine compound
    end

    def cosine compound
      RestClient.get(File.join(@@pug_proxy,"cid",@cid,"cosine",compound.cid)).to_f
    end
  end
end
