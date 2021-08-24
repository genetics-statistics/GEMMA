module TestHelpers

  # Runs gemma and returns true if successful
  def gemma(opts)
    system("./bin/gemma #{opts}")
  end

  def read(fn, line=0)
    count = 0
    File.open(fn, "r:utf-8").each_line { |ln|
      return ln.chomp.split("\t") if count == line
      count += 1
    }
  end

  def expect(fn, list)
    lines = File.read(fn).split("\n") # avoid this for large files
    list.each do | l |
      line,colnum,value = l
      cols = lines[line].chomp.split("\t")
      assert_equal value,cols[colnum]
    end
  end
end
