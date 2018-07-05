$1!="#" {s+=$3; s2+=$3*$3; field=$5; count++;}
END{
	print "# <s> \t error of <s> \t H \t count"
    	s = s/count; 
    	s2= sqrt(s2/count - s*s);

    	print s "\t" s2 "\t" field "\t" count
}
