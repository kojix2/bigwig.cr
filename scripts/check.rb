#!/usr/bin/env ruby
require "fileutils"
require "rbconfig"

Dir.chdir(File.expand_path("..", __dir__))
STDOUT.sync = true

def run(*command)
  puts "+ #{command.join(" ")}"
  return if system(*command)

  exit($?.exitstatus || 1)
end

run("crystal", "--version")
run("crystal", "spec", "--error-trace")
run("crystal", "run", "scripts/fixtures/generate.cr", "--", "--check")

FileUtils.mkdir_p("build")
Dir.glob("examples/*.cr").sort.each do |source|
  output = File.join("build", File.basename(source, ".cr"))
  run("crystal", "build", source, "--error-trace", "-o", output)
end

run("crystal", "run", "scripts/reference/generate.cr", "--", "build/interop")
run(RbConfig.ruby, "scripts/reference/verify.rb",
    "--interop-dir", "build/interop",
    "--report", "build/reference-validation.json")
