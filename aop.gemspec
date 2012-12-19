# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)

Gem::Specification.new do |s|
  s.name        = "aop"
  s.version     = File.read("./VERSION")
  s.authors     = ["Christoph Helma","Denis Gebele","Micha Rautenberg"]
  s.email       = ["helma@in-silico.ch","gebele@in-silico.ch","rautenenberg@in-silico.ch"]
  s.homepage    = "http://github.com/OpenTox/task"
  s.summary     = %q{AOP GUI}
  s.description = %q{AOP GUI}
  s.license     = 'GPL-3'
  #s.platform    = Gem::Platform::CURRENT

  s.rubyforge_project = "aop"

  s.files         = `git ls-files`.split("\n")
  #s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  #s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  #s.require_paths = ["lib"]
  s.required_ruby_version = '>= 1.9.2'

  # specify any dependencies here; for example:
  s.add_runtime_dependency "sinatra"
  s.add_runtime_dependency "sinatra-contrib"
  s.add_runtime_dependency "haml"
  s.add_runtime_dependency "yajl-ruby"
  s.add_runtime_dependency "rest-client"
  s.add_runtime_dependency "rest-client-components"
  s.add_runtime_dependency "memcache-client"
  s.add_runtime_dependency "unicorn"
end

