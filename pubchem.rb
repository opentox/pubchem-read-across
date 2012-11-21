require '../opentox-client/lib/opentox-client.rb'
require 'json'
require 'base64'

def Math.gauss(x, sigma = 0.3) 
  d = 1.0 - x.to_f
  Math.exp(-(d*d)/(2*sigma*sigma))
end

module OpenTox

  # doc @ http://pubchem.ncbi.nlm.nih.gov/pug_rest/
  class PubChemCompound < Compound
    attr_writer :cid
    attr_accessor :similarity, :p, :assays
    
    def initialize cid=nil
      @pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
      @cid = cid
      @assays = nil
      @similarity_threshold = 85
      @neighbors = nil
      @predicted_assays = nil
      #@predicted_targets = nil
      #@priors = {}
      #@priors = JSON.parse(File.read("priors.json"))
    end

    def fingerprint
      unless @fingerprint
        begin
          # ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
          base64key = `curl http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/#{cid}/SDF|grep -A1 PUBCHEM_CACTVS_SUBSKEYS|sed '1d'`.chomp
          @fingerprint = Base64.decode64(base64key)[4..-1].unpack("B*").first[0..-8].split(//).collect{|c| c == "1"}
        rescue
        end
      end
      @fingerprint
    end

    def self.from_name name
      pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
      compounds = []
      session[:name] = name
      cid = RestClientWrapper.get(File.join(pug_uri,URI.escape(name),"cids","TXT"))
      #puts response
      #response.split("\n") do |cid|
        puts cid
        compound = OpenTox::PubChemCompound.new
        compound.cid = cid.chomp
        compounds << compound
      #end
      compounds
    end

    def neighbors
      unless @neighbors
        @neighbors = []
        result = pubchem_search File.join(@pug_uri, "compound", "similarity", "cid", cid.to_s, "JSON")+"?Threshold=#{@similarity_threshold}&MaxRecords=100"
        while result["Waiting"] do
          sleep 2
          listkey = result["Waiting"]["ListKey"]
          result = pubchem_search File.join(@pug_uri, "compound", "listkey", listkey, "cids", "JSON")
          #result = pubchem_search File.join(@pug_uri, "compound", "listkey", listkey, "assaysummary", "JSON")
        end
        puts "Neighbor CIDs received"
        result["IdentifierList"]["CID"].each do |cid|
          unless cid.to_s == @cid.to_s
            c = PubChemCompound.new cid.to_s
            @neighbors << c if c.assays #and !(c.targets + c.non_targets).empty?
          end
        end if result and result["IdentifierList"]
=begin
        if result and result["Table"]
          columns = result["Table"]["Columns"]["Column"]
          table = result["Table"]["Row"].collect{|cell| cell.values.flatten}
          cid_idx = columns.index("CID")
          cids = table.collect{|r| r[cid_idx]}.uniq
          cids.each do |cid|
            unless cid.to_s == @cid.to_s
              tab = table.collect{|r| r if r[cid_idx] == cid}.compact
              c = PubChemCompound.new
              c.extract_result columns, tab
              c.similarity = tanimoto c
              @neighbors << c unless (c.targets + c.non_targets).empty?
            end
          end
        end
=end
        #@neighbors.sort!{|a,b| b.similarity <=> a.similarity}
      end
      @neighbors
    end

    def assays
      unless @assays
        result = pubchem_search File.join(@pug_uri, "compound", "cid", cid.to_s, "assaysummary", "JSON")
        extract_result result["Table"]["Columns"]["Column"], result["Table"]["Row"].collect{|cell| cell.values.flatten} if result and result["Table"]
      end
      @assays
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
      unless @predicted_assays
        @predicted_assays = []
        neighbors.collect{|n| n.assays.collect{|a| a["AID"]}}.flatten.compact.uniq.each do |aid|
          predicted_assay = {"AID" => aid}
          neighbors.each do |neighbor|
            if similarity(neighbor) > 0.5 # avoid downweighting
              search = neighbor.assays.select{|a| a["AID"] == aid}
              search.each do |assay|
                predicted_assay["Target GI"] ||= assay["Target GI"]
                predicted_assay["Target Name"] ||= assay["Target Name"]
                predicted_assay["Assay Name"] ||= assay["Assay Name"]
                predicted_assay[:active_similarities] ||= []
                predicted_assay[:inactive_similarities] ||= []

                if assay["Activity Outcome"] == "active"
                  predicted_assay[:p_active] ? predicted_assay[:p_active] = predicted_assay[:p_active]*similarity(neighbor) : predicted_assay[:p_active] = similarity(neighbor)
                  predicted_assay[:p_inactive] ? predicted_assay[:p_inactive] = predicted_assay[:p_inactive]*(1-similarity(neighbor)) : predicted_assay[:p_inactive] = 1-similarity(neighbor)
                  predicted_assay[:active_similarities] << similarity(neighbor)
                elsif assay["Activity Outcome"] == "inactive"
                  predicted_assay[:p_active] ? predicted_assay[:p_active] = predicted_assay[:p_active]*(1-similarity(neighbor)) : predicted_assay[:p_active] = 1-similarity(neighbor)
                  predicted_assay[:p_inactive] ? predicted_assay[:p_inactive] = predicted_assay[:p_inactive]*similarity(neighbor) : predicted_assay[:p_inactive] = similarity(neighbor)
                  predicted_assay[:inactive_similarities] << similarity(neighbor)
                end
              end
            end
          end
          if predicted_assay[:p_active] and predicted_assay[:p_inactive] and predicted_assay[:p_active] != 0 and predicted_assay[:p_inactive] != 0
            predicted_assay[:p_active] = predicted_assay[:p_active]/(predicted_assay[:p_active]+predicted_assay[:p_inactive])
            predicted_assay[:p_inactive] = predicted_assay[:p_inactive]/(predicted_assay[:p_active]+predicted_assay[:p_inactive])
            if predicted_assay[:p_active] > predicted_assay[:p_inactive]
              predicted_assay[:prediction] = "active"
            elsif predicted_assay[:p_active] < predicted_assay[:p_inactive]
              predicted_assay[:prediction] = "inactive"
            end
            @predicted_assays << predicted_assay
          end
        end
        #@predicted_targets.sort{|a,b| b[:p_active] <=> a[:p_active]}
      end
      @predicted_assays
    end

    def predicted_active_assays
      predicted_assays.select{|a| a[:prediction] == "active"} if predicted_assays
    end

    def predicted_inactive_assays
      predicted_assays.select{|a| a[:prediction] == "inactive"} if predicted_assays
    end

    def predicted_targets
      predicted_active_assays.select{|a| a[:target_gi]} if predicted_assays
    end

    def predicted_non_targets
      inactive_assays.select{|a| a[:target_gi]} if predicted_assays
    end

    def to_smiles
      RestClient.get(File.join(@pug_uri, "compound", "cid", cid.to_s, "property", "CanonicalSMILES", "TXT")).strip
    end

    def image_uri
      File.join @pug_uri, "compound", "cid", @cid, "PNG"#?record_type=3d&image_size=small"
    end

    def similarity compound
      cosine compound
    end

    def tanimoto compound
      if fingerprint and compound.fingerprint
        m11 = 0.0
        m1 = 0.0
        fingerprint.each_index do |i|
          m11 += 1 if (@fingerprint[i] and compound.fingerprint[i])
          m1 += 1 if (@fingerprint[i] or compound.fingerprint[i]) 
        end
        m11/m1
      end
    end

    def cosine compound
      if fingerprint and compound.fingerprint
        m11 = 0.0
        m01 = 0.0
        m10 = 0.0
        m00 = 0.0
        fingerprint.each_index do |i|
          m11 += 1 if (@fingerprint[i] and compound.fingerprint[i])
          m01 += 1 if (!@fingerprint[i] and compound.fingerprint[i]) 
          m10 += 1 if (@fingerprint[i] and !compound.fingerprint[i]) 
          m00 += 1 if (!@fingerprint[i] and !compound.fingerprint[i]) 
        end
        m11/((m01+m11)*(m10+m11))**0.5
      end
    end

=begin
      f1 = File.open(File.join(".","tmp",SecureRandom.uuid+".smi"),"w+")
      f1.puts to_smiles
      f1.close
      f2 = File.open(File.join(".","tmp",SecureRandom.uuid+".smi"),"w+")
      f2.puts compound.to_smiles
      f2.close
      sim = `babel #{f1.path} #{f2.path} -ofpt 2>/dev/null| grep Tanimoto|cut -d "=" -f2`.strip.to_f
      File.delete(f1.path)
      File.delete(f2.path)
      sim
    end
=end

    def pubchem_search url
      attempts = 0
      begin
        attempts += 1
        json = RestClient.get url, :timeout => 90000000
        puts url
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

    def extract_result columns, table
      @assays = []
      table.each do |row|
        @assays << {}
        row.each_with_index do |cell,i|
          if columns[i] == "CID" 
            @cid = cell if @cid.nil?
          else
            cell.blank? ?  @assays.last[columns[i]] = nil : @assays.last[columns[i]] = cell
          end
        end
      end
    end

    def priors aid
      unless @priors[aid]
        @priors[aid] = {"nr_active" => 0, "nr_inactive" => 0}
        result = nil
        result = pubchem_search File.join(@pug_uri, "assay", "aid", aid.to_s, "cids", "JSON?cids_type=active&list_return=listkey")
        @priors[aid]["nr_active"] = result["IdentifierList"]["Size"].to_i if result
        result = nil
        result = pubchem_search File.join(@pug_uri, "assay", "aid", aid.to_s, "cids", "JSON?cids_type=inactive&list_return=listkey")
        @priors[aid]["nr_inactive"] = result["IdentifierList"]["Size"].to_i if result
        File.open("priors.json","w+"){|f| f.puts @priors.to_json}
      end
      @priors[aid]
    end

=begin
    def assay_summary assay
        if assay["Target GI"] and !@assays[assay["AID"]]
          @assays[assay["AID"]] = {"nr_active" => 0, "nr_inactive" => 0}
          pubchem_search File.join(@pug_uri, "assay", "aid", assay["AID"].to_s, "cids", "JSON?cids_type=active")
          @assays[assay["AID"]]["nr_active"] = @result["InformationList"]["Information"].first["CID"].size if @result
          pubchem_search File.join(@pug_uri, "assay", "aid", assay["AID"].to_s, "cids", "JSON?cids_type=inactive")
          @assays[assay["AID"]]["nr_inactive"] = @result["InformationList"]["Information"].first["CID"].size if @result
          print "getting (in)actives for aid "
          puts assay["AID"]
          print @assays[assay["AID"]]["nr_active"]
          print " "
          puts @assays[assay["AID"]]["nr_inactive"]
          File.open("assays.json","w+"){|f| f.puts @assays.to_json}
        end
    end
=end

=begin

    def properties
      properties = [
        "XLogP",
        "ExactMass",
        "MonoisotopicMass",
        "TPSA",
        "Complexity",
        "Charge",
        "HBondDonorCount",
        "HBondAcceptorCount",
        "RotatableBondCount",
        "HeavyAtomCount",
        "IsotopeAtomCount",
        "AtomStereoCount",
        "DefinedAtomStereoCount",
        "UndefinedAtomStereoCount",
        "BondStereoCount",
        "DefinedBondStereoCount",
        "UndefinedBondStereoCount",
        "CovalentUnitCount",
        "Volume3D",
        "XStericQuadrupole3D",
        "YStericQuadrupole3D",
        "ZStericQuadrupole3D",
        "FeatureCount3D",
        "FeatureAcceptorCount3D",
        "FeatureDonorCount3D",
        "FeatureAnionCount3D",
        "FeatureCationCount3D",
        "FeatureRingCount3D",
        "FeatureHydrophobeCount3D",
        "ConformerModelRMSD3D",
        "EffectiveRotorCount3D",
        "ConformerCount3D",
      ]
      pubchem_search File.join(@pug_uri, "compound", "cid", @cid, "property", properties.join(","), "JSON")
      @result["PropertyTable"]["Properties"].first
    end

    def from_smiles smiles
      pubchem_search File.join(@pug_uri, "compound", "smiles", smiles, "assaysummary", "JSON")
      extract_result @result["Table"]["Columns"]["Column"], @result["Table"]["Row"].collect{|cell| cell.values.flatten}
    end
    def property_similarity compound
      svd = OpenTox::SVD.new(GSL::Matrix [[properties, compound.properties]])
      OpenTox::Algorithm::Similarity.cosine svd.data_transformed_matrix.first, svd.data_transformed_matrix.last
    end

    def assay_similarity compound
      tanimoto [[active_assays,inactive_assays],[compound.active_assays,compound.inactive_assays]]
    end

    def target_similarity compound
      tanimoto [[targets,non_targets],[compound.targets,compound.non_targets]]
    end

    def tanimoto features
      common = features.first.flatten & features.last.flatten
      same_outcome = (features.first.first & features.last.first) + (features.first.last & features.last.last)
      same_outcome.size.to_f/common.size
    end

    def euclid features
    end

    def to_name
      RestClient.get(File.join(@pug_uri, "compound", "cid", @cid, "property", "IUPACName", "TXT")).strip
    end
=end

  end

=begin
  class PubChemNeighbors < Dataset
    include PubChem

    attr_accessor :query, :neighbors

    def initialize
      @similarity_threshold = 95
      @neighbors = []
      @pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    end

    def from_smiles smiles
      #@query = PubChemCompound.new.from_smiles smiles
      pubchem_search File.join(@pug_uri, "compound", "similarity", "smiles", smiles, "JSON")+"?Threshold=#{@similarity_threshold}&MaxRecords=250"
      listkey = @result["Waiting"]["ListKey"]
      while @result["Waiting"] do
        sleep 1
        pubchem_search File.join(@pug_uri, "compound", "listkey", listkey, "assaysummary", "JSON")
      end
      #File.open("search.yaml","w+"){|s| s.puts @result.to_yaml}
      columns = @result["Table"]["Columns"]["Column"]
      table = @result["Table"]["Row"].collect{|cell| cell.values.flatten}
      cid_idx = columns.index("CID")
      cids = table.collect{|r| r[cid_idx]}.uniq
      cids.each do |cid|
        tab = table.collect{|r| r if r[cid_idx] == cid}.compact
        c = PubChemCompound.new
        c.extract_result columns, tab
        @neighbors << c unless (c.targets + c.active_assays).flatten.compact.empty?
      end
      @query = @neighbors.shift
      File.open("search.yaml","w+"){|s| s.puts self.to_yaml}
      #puts @neighbors.query.to_name
    end
  end
=end
end
