SERVICE = "aop"
require 'bundler'
Bundler.require
timeout = 600
require './application.rb'
run OpenTox::Application
