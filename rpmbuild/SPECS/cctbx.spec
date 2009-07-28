#-------------------------------------------------------------------------------
# Cameron Mura <cmura@virginia.edu>; 27jul09;
#-------------------------------------------------------------------------------
%define debug_package %{nil}
Name:     cctbx
Version:  2009_02_15_2320
Release:  3
Summary:  Computational Crystallography Toolbox
License:  Open Source (BSD-like)
Group:    Applications/Engineering
URL:	  http://cctbx.sourceforge.net
Source0:  %{name}_bundle_%{version}.tar.gz
Source1:  ccp4_mon_lib.tar.gz
Packager: Cameron Mura <cmura@virginia.edu>
#Patch0:  cctbx_bundle_StrucBio.patch
BuildRoot: /var/tmp/%{name}-buildroot

%description
The Computational Crystallography Toolbox (cctbx) is being developed as the
open source component of the PHENIX system. The goal of the PHENIX project is
to advance automation of macromolecular structure determination. PHENIX depends
on the cctbx, but not vice versa. This hierarchical approach enforces a clean
design as a reusable library. The cctbx is therefore also useful for
small-molecule crystallography and even general scientific applications. The
cctbx is designed with an open and flexible architecture to promote
extendability and easy incorporation into other software environments. The
package is organized as a set of ISO C++ classes with Python bindings. This
organization combines the computational efficiency of a strongly typed compiled
language with the convenience and flexibility of a dynamically typed scripting
language in a strikingly uniform and very maintainable way.  Note that CCTBX is
developed under an open source (BSD-like) license; for addition information,
see http://cctbx.svn.sourceforge.net/viewvc/cctbx/trunk/cctbx/LICENSE_2_0.txt.

%prep
%setup -a 0 -c

%build
./cctbx_install_script.csh /usr/bin/python 8
. ./cctbx_build/setpaths.sh
cp %{SOURCE1} .
libtbx.unpack_in_sources %{SOURCE1}

%install
[ $RPM_BUILD_ROOT != "/" ] && rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/usr/local/%{name}/%{version}
mkdir -p $RPM_BUILD_ROOT/etc/profile.d
echo '. "/usr/local/%{name}/%{version}/cctbx_build/setpaths.sh"' \
			 > $RPM_BUILD_ROOT/etc/profile.d/%{name}.sh
echo 'source "/usr/local/%{name}/%{version}/cctbx_build/setpaths.csh"' \
			 > $RPM_BUILD_ROOT/etc/profile.d/%{name}.csh
cp -R cctbx_sources cctbx_install_script.csh cctbx_build ${RPM_BUILD_ROOT}/usr/local/%{name}/%{version}/.
cp -R %{SOURCE1} ${RPM_BUILD_ROOT}/usr/local/%{name}/%{version}/.

%clean
[ $RPM_BUILD_ROOT != "/" ] && rm -rf $RPM_BUILD_ROOT

%post
# update paths and clean-up/regenerate byte-compiled (.pyc) python modules:
cd /usr/local/%{name}/%{version}/cctbx_build
/usr/bin/python /usr/local/%{name}/%{version}/cctbx_sources/libtbx/configure.py \
		fftw3tbx rstbx smtbx mmtbx clipper
cd /usr/local/%{name}/%{version}
find . -iname "*.pyc" -exec /bin/rm "{}" \;
/usr/local/cctbx/2009_02_15_2320/cctbx_build/bin/libtbx.py_compile_all 

%preun

%files
%defattr(-,root,root)
#%doc ChangeLog.txt ChangeLog.html INSTALL* NOTICE README*
/usr/local/%{name}/%{version}
/etc/profile.d/%{name}.sh
/etc/profile.d/%{name}.csh

%changelog
* Tue Jul 28 2009 Cameron Mura <cmura@virginia.edu) 2009_02_15_2320-3
- Very minor clean-up before checking into cctbx's svn repo...

* Mon Jul 27 2009 Cameron Mura <cmura@virginia.edu) 2009_02_15_2320-2
- Revise the build procedure to make use of cctbx_sources/libtbx/configure.py,
  rather than the ugly (and potentially dangerous) approach of building in the
  system-wide /usr/local/ space.

* Sun Jul 26 2009 Cameron Mura <cmura@virginia.edu) 2009_02_15_2320-1
- Initial RPM build on Fedora 10.

