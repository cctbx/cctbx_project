def touch_init(env, target, source):
  assert len(target) == 1
  print "Creating:", str(target[0])
  open(str(target[0]), "w").close()
