require 'sinatra'
require "sinatra/reloader" 
require "haml"
require 'yaml'
require "./pubchem.rb"
also_reload './pubchem.rb'

get '/?' do
  #@neighbors = OpenTox::PubChemNeighbors.new
  smiles = "OC(=O)C1=C(C=CC=C1)OC(=O)C"
  #smiles = "c1cc(CC)ccc1"
  #smiles = "CC(=O)Nc1ccc(O)cc1"
  smiles = "C1=CC(=C(C=C1Cl)Cl)OCC(=O)O"
  #@compound.from_smiles smiles
  #@neighbors.from_smiles smiles
=begin
  File.open("compound.yaml","w+"){|f| f.puts @compound.to_yaml}
  @compound = YAML.load_file "compound.yaml"
=end
  @neighbors = YAML.load_file "search.yaml"
  haml :index
end
