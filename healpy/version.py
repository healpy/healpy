rev  = "$Revision$"

revnumber = int(rev.strip('$')[9:])
date = "$Date$"

__version__='trunk-r%d'%(revnumber)
