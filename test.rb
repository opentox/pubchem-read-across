require 'test/unit'
#require 'yaml'
require './pubchem.rb'

class MRATest < Test::Unit::TestCase

  def setup
    #@p = PubChem::Compound.new "c1cc(CC)ccc1"
    @p = PubChem::Compound.new
    @p.from_smiles "CC(=O)Nc1ccc(O)cc1"
  end

  def test_initialize
    puts @p.active_assays.size
    puts @p.inactive_assays.size
    puts @p.targets.size
    puts @p.non_targets.size
    assert_equal "7500", @p.cid
    #assert_equal true, @p.aids[:active].include?(1188)
    #assert_equal true, @p.aids[:inactive].include?(435)
  end

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

end
