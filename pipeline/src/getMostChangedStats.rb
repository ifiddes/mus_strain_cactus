#!/usr/bin/env ruby
require 'nokogiri'

(warn "Usage: #{$0} statsXML1 statsXML2"; exit(1)) if ARGV.size != 2

MIN_TRANSCRIPT_FRAC=0.02
def get_transcripts_per_tag(xml)
  h = {}
  xml.traverse { |node| h[node.path] = node["transcripts"].to_i }
  h
end

def normalize(transcripts_per_tag)
  ret = {}
  total_transcripts = transcripts_per_tag["/stats/ok"] + transcripts_per_tag["/stats/not_ok"]
  transcripts_per_tag.each {|k,v| ret[k] = v.to_f/total_transcripts }
  ret
end

stats1 = Nokogiri::XML(open(ARGV[0]))
stats2 = Nokogiri::XML(open(ARGV[1]))

transcripts_per_tag_1 = get_transcripts_per_tag(stats1)
transcripts_per_tag_2 = get_transcripts_per_tag(stats2)

transcript_frac_per_tag_1 = normalize(transcripts_per_tag_1).select { |k, v| v >= MIN_TRANSCRIPT_FRAC }
transcript_frac_per_tag_2 = normalize(transcripts_per_tag_2).select { |k, v| v >= MIN_TRANSCRIPT_FRAC }

puts "tag\tfoldChange\ttranscripts1\ttranscripts2\ttranscriptFrac1\ttranscriptFrac2"
foldChanges = {}
transcript_frac_per_tag_1.each { |k, v| foldChanges[k] = transcript_frac_per_tag_2[k] ? transcript_frac_per_tag_2[k]/v : 0 }
transcript_frac_per_tag_2.each { |k, v| foldChanges[k] = Float::INFINITY unless foldChanges[k] }
foldChanges.to_a.sort {|e1, e2| e1[1] <=> e2[1]}.each do |e|
  tag, foldChange = e
  puts "#{tag}\t#{foldChange}\t#{transcripts_per_tag_1[tag]}\t#{transcripts_per_tag_2[tag]}\t#{transcript_frac_per_tag_1[tag]}\t#{transcript_frac_per_tag_2[tag]}"
end
