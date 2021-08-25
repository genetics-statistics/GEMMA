module TestHelpers

  # Runs gemma and returns true if successful
  def gemma(opts)
    assert system("./bin/gemma #{opts}")
  end

  def read(fn, line=0)
    count = 0
    File.open(fn, "r:utf-8").each_line { |ln|
      return ln.chomp.split("\t") if count == line
      count += 1
    }
  end

  def expect(fn, list)
    lines = File.read(fn).split("\n")
    lines = lines.map { |l| l.split("\t") } # avoid this for large files
    list.each do | l |
      line,colnum,value = l
      if colnum.is_a? String
        colnum = lines[0].index(colnum)
      end
      cols =
        if line == :max
          lines.max_by {|a| a[colnum].to_f}
        else
          lines[line]
        end
      # assert_equal value,cols[colnum]
      assert_in_delta value.to_f,cols[colnum].to_f, 0.001
    end
  end
end
