#!/usr/bin/env ruby
require "./pubchem.rb"

false_positives = YAML.load_file("false_positives.yaml")
#false_positives = []
until false_positives.size > 100 do
  result = {}
  @compound = OpenTox::PubChemCompound.new
  # http://www.ncbi.nlm.nih.gov/sites/entrez?term=all%5Bfilt%5D&cmd=search&db=pccompound
  @compound.cid = Random.new.rand(1..35611104)
  puts @compound.cid
  if @compound.targets and @compound.non_targets and !(@compound.targets + @compound.non_targets).empty?
    puts "predicting ..."
    result[:cid] = @compound.cid
    measured_non_targets = @compound.non_targets.collect{|t| t["Target GI"]}.compact.uniq
    predicted_targets = @compound.predicted_targets.collect{|t| t[:target_gi] if t[:prediction] == "active"}.compact.uniq

    result[:fp_targets] = {}
    (predicted_targets & measured_non_targets).each do |gi|
      result[:fp_targets][gi] = {:p => {:active => nil, :inactive => nil}, :measured => [], :neighbors => []}
      result[:fp_targets][gi][:measured] = @compound.inactive_assays.select{|a| a["Target GI"] == gi}
      result[:fp_targets][gi][:p][:active] = @compound.predicted_targets.collect{|t| t[:p_active] if t[:target_gi] == gi}.compact.uniq
      result[:fp_targets][gi][:p][:inactive] = @compound.predicted_targets.collect{|t| t[:p_inactive] if t[:target_gi] == gi}.compact.uniq 

      @compound.neighbors.select{|n| n.assays.collect{|a| a["Target GI"]}.include? gi }.each do |neighbor|
        result[:fp_targets][gi][:neighbors] << {
          :cid => neighbor.cid,
          :similarity => neighbor.similarity,
          :assays => neighbor.assays.select{|a| a["Target GI"] == gi }
        }
      end
    end
    unless result[:fp_targets].empty?
      false_positives << result
      File.open("false_positives.yaml","w+") {|f| f.puts false_positives.to_yaml}
    end
  end
end

