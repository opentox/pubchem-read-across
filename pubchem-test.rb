require 'test/unit'
require '../opentox-client/lib/opentox-client'
require './pubchem.rb'
require 'yaml'

class AOPTest < Test::Unit::TestCase

  def setup
    @pug_uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    @compound = OpenTox::PubChemCompound.new #3036
    @compound.cid = 1983
    #@compound.from_name "2,4-D"
  end

  def test_initialize
    print @compound.cid
    print " "
    puts @compound.to_smiles
    puts "measured targets"
    puts @compound.targets.collect{|t| t["Target Name"]}.to_yaml
=begin
    puts "predicted targets"
    puts @compound.predicted_targets.select{|t| t[:prediction] == "active"}.size

    puts "predicted non_targets"
    puts @compound.predicted_targets.select{|t| t[:prediction] == "inactive"}.size
    #puts @compound.predicted_non_targets.values.inspect
    measured_target_gis = @compound.targets.collect{|t| t["Target GI"]}.compact.uniq
    measured_nontarget_gis = @compound.non_targets.collect{|t| t["Target GI"]}.compact.uniq
    predicted_target_gis = @compound.predicted_targets.collect{|t| t[:target_gi] if t[:prediction] == "active"}.compact.uniq
    predicted_nontarget_gis = @compound.predicted_targets.collect{|t| t[:target_gi] if t[:prediction] == "inactive"}.compact.uniq
    print "correct predicted targets: "
    puts (measured_target_gis & predicted_target_gis).size
    print "new predicted targets: "
    puts (predicted_target_gis - measured_target_gis).size
    print "correct predicted non-targets: "
    puts (measured_nontarget_gis & predicted_nontarget_gis).size
    print "new predicted non-targets: "
    puts (predicted_nontarget_gis - measured_nontarget_gis).size
    print "incorrect predicted targets: "
    puts (measured_nontarget_gis & predicted_target_gis).size
    puts (measured_nontarget_gis & predicted_target_gis).sort.to_yaml
    puts @compound.predicted_targets.select{|t| t[:prediction] == "active"}.to_yaml
    print "incorrect predicted non-targets: "
    puts (measured_target_gis & predicted_nontarget_gis).size
=end
=begin
    @compound.neighbors.each do |n|
      #print n.cid
      #print " "
      print n.to_smiles
      print " "
      print n.similarity
      print " "
      puts n.p
      puts n.targets.sort.inspect
      #puts n.non_targets.inspect
    end
=end
    #File.open("Acetaminophen.yaml","w+"){|f| f.puts @compound.to_yaml}
    #puts @compound.neighbors.size
    #assert_equal "7500", @p.cid
    #assert_equal true, @p.aids[:active].include?(1188)
    #assert_equal true, @p.aids[:inactive].include?(435)
  end
=begin
  def test_similarity_search
    #puts @p.to_smiles
    @p.neighbors.each do |n|
      #puts n.to_smiles
      puts @p.target_similarity(n)
    end
    #puts @p.neighbors.inspect
    #assert_equal 100, @p.neighbor_cids.size
  end

  def test_assay_description
    puts @p.assay_description.to_yaml
  end

  def test_assay_genes
    puts @p.assay_genes.to_yaml
  end

  def test_assay_similarity
    @p2 = PubChem::Compound.new "OC(=O)C1=C(C=CC=C1)OC(=O)C"
    puts @p.assay_similarity(@p2)
  end

  def test_target_similarity
    @p2 = PubChem::Compound.new "OC(=O)C1=C(C=CC=C1)OC(=O)C"
    puts @p.target_similarity(@p2)
  end
=end

end
