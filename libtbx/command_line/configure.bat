@setlocal
@for /F "delims=" %%i in ('libtbx.show_dist_paths libtbx') do @set d=%%i
@"%LIBTBX_PYEXE%" "%d%\env_config.py" %*
