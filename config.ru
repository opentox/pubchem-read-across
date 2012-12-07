SERVICE = "aop"
require 'bundler'
Bundler.require
require './application.rb'
use Rack::Session::Dalli, :cache => Dalli::Client.new
#use Rack::Cache,
#  :verbose => true,
#  :metastore   => "memcached://127.0.0.1:11211/meta",
#  :entitystore => "memcached://127.0.0.1:11211/body"
run OpenTox::Application
