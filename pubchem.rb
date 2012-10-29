require '../opentox-client/lib/opentox-client.rb'
require 'json'

def Math.gauss(x, sigma = 0.3) 
  d = 1.0 - x.to_f
  Math.exp(-(d*d)/(2*sigma*sigma))
end

module PubChem

  attr_accessor :result

  def initialize
    @pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
  end

  def pubchem_search url
    json =  RestClient.get url#, :accept => "application/json", :timeout => 90000000
    @result = JSON.parse json
  rescue
    puts url
    puts $!.message
    @result = nil
  end

end

module OpenTox

  class PubChemCompound < Compound
    include PubChem
    # doc @ http://pubchem.ncbi.nlm.nih.gov/pug_rest/
    attr_writer :cid
    attr_accessor :similarity, :p, :assays
    
    def initialize 
      super
      @summary = []
      @similarity_threshold = 75
      @neighbors = []
      @predicted_targets = []
    end

    def from_name name
      @inchi = RestClientWrapper.get File.join(CACTUS_URI,URI.escape(name),"stdinchi")
    end

    def neighbors
      if @neighbors.empty?
        pubchem_search File.join(@pug_uri, "compound", "similarity", "cid", cid.to_s, "JSON")+"?Threshold=#{@similarity_threshold}&MaxRecords=100"
        listkey = @result["Waiting"]["ListKey"]
        while @result["Waiting"] do
          sleep 1
          pubchem_search File.join(@pug_uri, "compound", "listkey", listkey, "assaysummary", "JSON")
        end
        columns = @result["Table"]["Columns"]["Column"]
        table = @result["Table"]["Row"].collect{|cell| cell.values.flatten}
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
        @neighbors.sort!{|a,b| b.similarity <=> a.similarity}
      end
      @neighbors
    end

    def summary
      if @summary.empty?
        pubchem_search File.join(@pug_uri, "compound", "cid", cid.to_s, "assaysummary", "JSON")
        extract_result @result["Table"]["Columns"]["Column"], @result["Table"]["Row"].collect{|cell| cell.values.flatten}
      end
      @summary
    end

    def active_assays
      summary.select{|a| a["Activity Outcome"] == "active"}
    end

    def inactive_assays
      summary.select{|a| a["Activity Outcome"] == "inactive"}
    end

    def targets
      active_assays.select{|a| a["Target GI"]}
    end

    def non_targets
      inactive_assays.select{|a| a["Target GI"]}
    end

    def predicted_targets
      if @predicted_targets.empty?
        target_gis = neighbors.collect{|n| n.summary.collect{|a| a["Target GI"]}}.flatten.compact.uniq
        target_gis.each do |gid|
          target = {:target_gi => gid}
          neighbors.each do |neighbor|
            if neighbor.similarity > 0.5 # avoid downweighting
              search = neighbor.summary.select{|a| a["Target GI"] == gid}
              unless search.empty? or search.size == 1
                print "+++ ("
                print search.size
                puts ")"
                puts search.inspect
              end
              search.each do |assay|
                target[:aid] ||= assay["AID"]
                target[:name] ||= assay["Target Name"]
                target[:assay_name] ||= assay["Assay Name"]
                target[:active_similarities] ||= []
                target[:inactive_similarities] ||= []

                if assay["Activity Outcome"] == "active"
                  target[:p_active] ? target[:p_active] = target[:p_active]*neighbor.similarity : target[:p_active] = neighbor.similarity
                  target[:p_inactive] ? target[:p_inactive] = target[:p_inactive]*(1-neighbor.similarity) : target[:p_inactive] = 1-neighbor.similarity
                  target[:active_similarities] << neighbor.similarity
                elsif assay["Activity Outcome"] == "inactive"
                  target[:p_active] ? target[:p_active] = target[:p_active]*(1-neighbor.similarity) : target[:p_active] = 1-neighbor.similarity
                  target[:p_inactive] ? target[:p_inactive] = target[:p_inactive]*neighbor.similarity : target[:p_inactive] = neighbor.similarity
                  target[:inactive_similarities] << neighbor.similarity
                end
              end
            end
          end
          if target[:p_active] and target[:p_inactive] and target[:p_active] + target[:p_inactive] != 0
            target[:p_active] = target[:p_active]/(target[:p_active]+target[:p_inactive])
            target[:p_inactive] = target[:p_inactive]/(target[:p_active]+target[:p_inactive])
            if target[:p_active] > target[:p_inactive]
              target[:prediction] = "active"
            elsif target[:p_active] < target[:p_inactive]
              target[:prediction] = "inactive"
            end
            @predicted_targets << target
          end
        end
        @predicted_targets.sort{|a,b| b[:p_active] <=> a[:p_active]}
      end
      @predicted_targets
    end

    def to_smiles
      RestClient.get(File.join(@pug_uri, "compound", "cid", cid.to_s, "property", "CanonicalSMILES", "TXT")).strip
    end

    def tanimoto compound
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

    def extract_result columns, table
      table.each do |row|
        @summary << {}
        row.each_with_index do |cell,i|
          if columns[i] == "CID" 
            @cid = cell if @cid.nil?
          else
            cell.blank? ?  @summary.last[columns[i]] = nil : @summary.last[columns[i]] = cell
          end
        end
      end
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

    def to_image_uri
      File.join @pug_uri, "compound", "cid", @cid, "PNG?record_type=3d&image_size=small"
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
