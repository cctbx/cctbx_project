from cctbx.web.asu_gallery import web_links
from cctbx.web.asu_gallery import html_head_title
import sys

def write_html(f=None):
  if (f is None): f = sys.stdout
  iucrcompcomm_jul2003 = web_links.iucrcompcomm_jul2003
  title = "ASU Gallery - Guide to notation"
  print >> f, html_head_title(title=title)
  print >> f, """\
<body>

<hr>
<h1>%(title)s</h1>
<hr>
Reference:
<a href="%(iucrcompcomm_jul2003)s"
>IUCr Computing Commission Newsletter No. 2, July 2003</a>
<hr>

Each cut plane of an asymmetric unit is defined by a condition of
the form
<pre>
  h*x+k*y+l*z+c>=0
</pre>
<tt>h</tt>,<tt>k</tt>,<tt>l</tt> are Miller indices and define the
normal vector of the cut plane, <tt>c</tt> is a constant which determines
the distance from the origin. <tt>x</tt>,<tt>y</tt>,<tt>z</tt> are
fractional coordinates in direct space. The expression
<tt>h*x+k*y+l*z+c</tt> is

<ul>
<li>exactly zero for points in the cut plane.
<li>greater than zero for points inside the asymmetric unit.
<li>less than zero for points outside the asymmetric unit.
</ul>

If all points that are exactly in a cut plane are not inside the
asymmetric unit, the condition changes from <tt>h*x+k*y+l*z+c&gt;=0</tt>
to <tt>h*x+k*y+l*z+c&gt;0</tt>.

<p>

To enhance readability the asymmetric unit conditions (shown under the pictures
in the gallery) are simplified by omitting terms with zeros (e.g.
<tt>0*x</tt>) and unit factors (e.g. <tt>x</tt> instead of
<tt>1*x</tt>). The constant term <tt>c</tt> is moved to the right-hand
side. For example:

<pre>
  x>=0
  y<=1/4
  z<1
  x-y<=1/2
</pre>

A point <tt>x</tt>,<tt>y</tt>,<tt>z</tt> is inside the asymmetric unit
only if all conditions are simultaneously true.

<p>

Often a face or edge on the surface of the asymmetric unit is only partially
inside. The dividing lines are defined by face- or edge-specific
sub-conditions. For example:

<pre>
  y<=1/4 [z<=1/2]
</pre>

The first condition defines the face as before. The second condition
in square brackets only applies if <tt>y=1/4</tt>.
This notation is recursive. For example:

<pre>
  y<=1/4 [z<=1/2 [x<=1/4]]
</pre>

The third condition only applies if <tt>y=1/4</tt> and <tt>z=1/2</tt>.

<p>
Some asymmetric units require the combination of conditions with
the boolean operators <i>and</i> or <i>or</i>. For example:

<pre>
  y>=0 [x<=0 | x>=1/4]
  y<=1/4 [z>=1/8 & z<=5/8]
</pre>

In words:

<ul>
<li>If <tt>y</tt> is exactly zero, a point is inside the asymmetric
    unit only if <tt>x</tt> is less than or equal to zero or
    greater than or equal to <tt>1/4</tt>.
<li>If <tt>y</tt> is exactly <tt>1/4</tt>, a point is inside the
    asymmetric unit only if <tt>z</tt> is greater than or equal to
    <tt>1/8</tt> and at the same time less than or equal to <tt>5/8</tt>.
</ul>

The boolean operators may occur at any level in the recursive
hierarchy defined by nested square brackets. An example is the
asymmetric unit of <a href="asu_088.html">space group
I&nbsp;41/a&nbsp;(No.&nbsp;88)</a>.

<hr>
<a href="index.html">Gallery of direct-space asymmetric units</a>

<hr>
</body>
</html>""" % vars()

if (__name__ == "__main__"):
  write_html()
