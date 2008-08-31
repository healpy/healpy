/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Class for storing parameter information for later use
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#include "fitshandle.h"
#include "simparams.h"

using namespace std;

void simparams::add_keys (ostream &os) const
  {
  for (unsigned int m=0; m<source_files.size(); ++m)
    {
    os << "ancestor"+dataToString(m+1)+"="+source_files[m] << endl;
    }
  vector<Param>::const_iterator iter = paramMap.begin();
  while (iter!=paramMap.end())
    {
    if (iter->comment != "")
      os << "# "+iter->comment << endl;
    if (iter->key != "")
      os << iter->key << "=" << iter->value << endl;
    ++iter;
    }
  }

void simparams::add_keys (fitshandle &out) const
  {
  fitshandle inp;
  for (unsigned int m=0; m<source_files.size(); ++m)
    {
    inp.open(source_files[m]);
    inp.goto_hdu(hdus[m]);
    out.add_comment("imported from HDU "+dataToString(hdus[m])+" of");
    out.add_comment(source_files[m]);
    out.copy_header(inp);
    out.add_comment("End of imported HDU");
    inp.close();
    }

  vector<Param>::const_iterator iter = paramMap.begin();
  while (iter!=paramMap.end())
    {
    if (iter->shortkey != "")
      out.add_key (iter->shortkey, iter->value, iter->comment);
    else
      out.add_comment (iter->comment);
    ++iter;
    }
  }
