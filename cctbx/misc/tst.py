from cctbx_boost import dev
def run():
  dev.waiting_loop(500, 1000000)
if (__name__ == "__main__"):
  import profile
  profile.run("run()")
