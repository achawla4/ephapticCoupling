function capacity = channelcapacityBSC(p)
Hbp = -p*log2(p)-(1-p)*log2(1-p);
capacity = 1 - Hbp;