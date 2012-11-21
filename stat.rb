#!/usr/bin/env ruby
require 'yaml'

stat = {:tp => 0, :tn => 0, :fp => 0, :fn => 0, :tp_p => [], :fp_p => []}
thresh = 0.05
Dir["./validation/*yaml"].each do |f|
  data = YAML.load_file f
  pa = data[:predicted][:active].select{|gi| data[:predicted][:p][gi][:p_active] > thresh }
  pi = data[:predicted][:inactive].select{|gi| data[:predicted][:p][gi][:p_inactive] > thresh }
  stat[:tp] += (pa & data[:measured][:active]).size
  stat[:tn] += (pi & data[:measured][:inactive]).size
  stat[:fp] += (pa & data[:measured][:inactive]).size
  stat[:fn] += (pi & data[:measured][:active]).size
  (pa & data[:measured][:active]).each{|gi| stat[:tp_p] << data[:predicted][:p][gi][:p_active] }
  (pa & data[:measured][:inactive]).each{|gi| stat[:fp_p] << data[:predicted][:p][gi][:p_active] }
end

stat[:tp_p].sort!
stat[:fp_p].sort!
puts stat.to_yaml
print "accuracy: "
puts (stat[:tp]+stat[:tn])/(stat[:tp]+stat[:tn]+stat[:fp]+stat[:fn]).to_f
print "sensitivity: "
puts stat[:tp]/(stat[:tp]+stat[:fn]).to_f
print "specificity: "
puts stat[:tn]/(stat[:tn]+stat[:fp]).to_f
print "positive predictive value: "
puts stat[:tp]/(stat[:tp]+stat[:fp]).to_f
print "negative predictive value: "
puts stat[:tn]/(stat[:tn]+stat[:fn]).to_f
