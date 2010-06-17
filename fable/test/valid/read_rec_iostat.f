      program prog
      dimension nums(2)
      read(11, rec=21, iostat=ios) num
      read(12, iostat=ios, err=10, end=20) num
      read(13, rec=23, iostat=ios) (nums(i), i=1,2)
      read(14, iostat=ios, err=30, end=40) (nums(i), i=1,2)
   10 continue
   20 continue
   30 continue
   40 continue
      end
