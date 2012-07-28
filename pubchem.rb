require '../opentox-client/lib/opentox-client.rb'
require 'json'

# get assay from endpoint
# search in endpoint ontology
#
# get measurements
# search for compound and assay 
#
# get affected pathways
# search for compound and genes
# identify affected pathways
# identify relations between affected genes/pathways and endpoint
#
# get related assays
# search for assays in ontology tree
# search for compound and related assays
#
# get similar compounds
# search for similar compounds

module PubChem

  def pubchem_search url
    puts url
    #json =  RestClient.get url, :accept => "application/json", :timeout => 90000000
    json =  `curl "#{url}"`#, :accept => "application/json", :timeout => 90000000
    @result = JSON.parse json
  end

  class Assay
    attr_accessor :aid
  end

  class Result
    attr_accessor :aid, :cid, :sid
  end

  class Substance
  end

  class Compound
    # doc @ http://pubchem.ncbi.nlm.nih.gov/pug_rest/
    include OpenTox
    include PubChem
    attr_accessor :result, :cid, :neighbors, :tanimoto
    
    def initialize 
      @uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
      @similarity_threshold = 95
      @summary = []
      @neighbors = []
    end

    def to_smiles
      RestClient.get(File.join(@uri, "compound", "cid", @cid, "property", "CanonicalSMILES", "TXT")).strip
    end

    def name
      RestClient.get(File.join(@uri, "compound", "cid", @cid, "property", "IUPACName", "TXT")).strip
    end

    def image
      File.join @uri, "compound", "cid", @cid, "PNG?record_type=3d&image_size=small"
    end

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
      pubchem_search File.join(@uri, "compound", "cid", @cid, "property", properties.join(","), "JSON")
      @result["PropertyTable"]["Properties"].first
    end

    def from_smiles smiles
      pubchem_search File.join(@uri, "compound", "smiles", smiles, "assaysummary", "JSON")
      from_summary @result["Table"]["Columns"]["Column"], @result["Table"]["Row"].collect{|cell| cell.values.flatten}
    end

    def from_summary columns, table
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

    def active_assays
      @summary.collect{|a| a if a["Activity Outcome"] == "active"}.compact
    end

    def inactive_assays
      @summary.collect{|a| a if a["Activity Outcome"] == "inactive"}.compact
    end

    def targets
      active_assays.collect{|a| a["Target Name"]}.compact
    end

    def non_targets
      inactive_assays.collect{|a| a["Target Name"]}.compact
    end

    def assay_similarity compound
      a1 = active_assays.collect{|a| a["Assay Name"]}
      a2 = compound.active_assays.collect{|a| a["Assay Name"]}
      i1 = inactive_assays.collect{|a| a["Assay Name"]}
      i2 = compound.inactive_assays.collect{|a| a["Assay Name"]}
      self_assays = a1 + i1
      compound_assays = a2 + i2
      common_assays = self_assays & compound_assays
      same_outcome = (a1 & a2) + (i1 & i2)
      same_outcome.size.to_f/common_assays.size
    end

    def target_similarity compound
      self_assays = targets + non_targets
      compound_assays = compound.targets + compound.non_targets
      common_assays = self_assays & compound_assays
      same_outcome = (targets & compound.targets) + (non_targets & compound.non_targets)
      same_outcome.size.to_f/common_assays.size
    end

    def assay_genes
      active = []
      @aids[:active].each do |aid|
        begin
        pubchem_search File.join(@uri, "assay", "aid", aid.to_s, "genes", "JSON")
        active << @result["InformationList"]["Information"].collect{|i| i["GeneID"]}.flatten
        rescue; end
      end
      @aids[:inactive].each do |aid|
        begin
        pubchem_search File.join(@uri, "assay", "aid", aid.to_s, "genes", "JSON")
        inactive << @result["InformationList"]["Information"].collect{|i| i["GeneID"]}.flatten
        rescue; end
      end
      {:active => active, :inactive => inactive } 
    end

    def get_neighbors smiles
      pubchem_search File.join(@uri, "compound", "similarity", "smiles", smiles, "JSON")+"?Threshold=#{@similarity_threshold}&MaxRecords=250"
      listkey = @result["Waiting"]["ListKey"]
      while @result["Waiting"] do
        sleep 1
        pubchem_search File.join(@uri, "compound", "listkey", listkey, "assaysummary", "JSON")
      end
      File.open("search.yaml","w+"){|s| s.puts @result.to_yaml}
      columns = @result["Table"]["Columns"]["Column"]
      table = @result["Table"]["Row"].collect{|cell| cell.values.flatten}
      cid_idx = columns.index("CID")
      cids = table.collect{|r| r[cid_idx]}.uniq
      cids.each do |cid|
        tab = table.collect{|r| r if r[cid_idx] == cid}.compact
        c = PubChem::Compound.new
        c.from_summary columns, tab
        @neighbors << c unless (c.targets + c.active_assays).flatten.compact.empty?
      end
      File.open("smiles.smi","w+"){|f| f.puts @neighbors.collect{|n| n.to_smiles}.join("\n")}
      `babel smiles.smi -ofpt 2>/dev/null| grep Tanimoto|cut -d "=" -f2`.split("\n").each_with_index do |t,i|
        @neighbors[i].tanimoto = t.strip.to_f
      end
    end
  end
end

