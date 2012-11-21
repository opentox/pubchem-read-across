#!/usr/bin/env ruby
require "./pubchem.rb"

until Dir["./validation/*.yaml"].size > 1000 do
#100.times do
  result = {}
  @compound = OpenTox::PubChemCompound.new
  # http://www.ncbi.nlm.nih.gov/sites/entrez?term=all%5Bfilt%5D&cmd=search&db=pccompound
  @compound.cid = Random.new.rand(1..35611104)
  puts @compound.cid
  unless File.exists?  "./validation/#{@compound.cid}.yaml"
    if (@compound.targets + @compound.non_targets).size > 0
    #if @compound.assays#.empty?
      begin
        puts "predicting ..."
        result[:cid] = @compound.cid
        result[:measured] = {}
        result[:predicted] = {}
        result[:measured][:active] = @compound.targets.collect{|t| t["Target GI"]}.compact.uniq
        result[:measured][:inactive] = @compound.non_targets.collect{|t| t["Target GI"]}.compact.uniq
        result[:predicted][:active] = @compound.predicted_targets.collect{|t| t[:target_gi] if t[:prediction] == "active"}.compact.uniq
        result[:predicted][:inactive] = @compound.predicted_targets.collect{|t| t[:target_gi] if t[:prediction] == "inactive"}.compact.uniq
        result[:predicted][:p] = {}
        @compound.predicted_targets.each do |t|
          result[:predicted][:p][t[:target_gi]] = {}
          result[:predicted][:p][t[:target_gi]][:p_active] = t[:p_active]
          result[:predicted][:p][t[:target_gi]][:p_inactive] = t[:p_inactive]
        end
        File.open("./validation/#{@compound.cid}.yaml","w+"){|f| f.puts result.to_yaml}
        puts result.to_yaml
      rescue
      end
    end
  end
end
