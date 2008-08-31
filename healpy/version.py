rev  = "$Revision: 112 $"

revnumber = int(rev.strip('$')[9:])
date = "$Date: 2008-07-20 17:46:55 +0200 (dim, 20 jui 2008) $"

__version__='trunk-r%d'%(revnumber)
