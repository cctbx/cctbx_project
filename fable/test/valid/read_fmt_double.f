      program prog
      double precision val
      character buf*9
      buf = '1234'
      read(buf, fmt='(bnf3.0)') val
      write(6, *) val
      read(buf, fmt='(bnf3.1)') val
      write(6, *) val
      read(buf, fmt='(bnf3.2)') val
      write(6, *) val
      buf = '.234'
      read(buf, fmt='(bnf3.2)') val
      write(6, *) val
      buf = '1.34'
      read(buf, fmt='(bnf3.2)') val
      write(6, *) val
      buf = '12.4'
      read(buf, fmt='(bnf3.2)') val
      write(6, *) val
      buf = '1 34'
      read(buf, fmt='(bnf3.0)') val
      write(6, *) val
      read(buf, fmt='(bzf3.0)') val
      write(6, *) val
      buf = '1234'
      read(buf, fmt='(1pbnf3.0)') val
      write(6, *) val
      read(buf, fmt='(2pbnf3.0)') val
      write(6, *) val
      read(buf, fmt='(-1pbnf3.0)') val
      write(6, *) val
      read(buf, fmt='(-2pbnf3.0)') val
      write(6, *) val
      buf = '1e34'
      read(buf, fmt='(1pbnf3.0)') val
      write(6, *) val
      buf = '3 5 e 1 2'
      read(buf, fmt='(bnf9.0)') val
      write(6, *) val
      read(buf, fmt='(bzf9.0)') val
      write(6, *) val
      read(buf, fmt='(bnf9.1)') val
      write(6, *) val
      read(buf, fmt='(bzf9.1)') val
      write(6, *) val
      end
