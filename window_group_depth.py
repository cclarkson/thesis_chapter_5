sync=open('/data/gambiae/chris3/Arabiensis_GWAS/Igor_aligned/Arabiensis_GWAS_with_merus/All_Igor_with_merus_raw_sync/Arab_Igor_with_Merus_2R.sync Arab_Igor_with_Merus_2R.sync', 'r')
out=open('MoshiA_MoshiD_Tarime_depth_2R_HEAD', 'w')
#out = open('Merus_depth', 'w')


#N=2
#for i in range(N): #check head of file to make sure correct
 #   line=sync.readline().strip()
 #   print(line)

out.write("chr\tposition\tMoshiAlive\tMoshiDead\tTarime\n")

for position in sync:
    position = position.rstrip()
    lineSplit = position.split()
    chrom = lineSplit[0]
    where = lineSplit[1]
    justAlleles=lineSplit[3:]
    moshiAlive = justAlleles[0:10]
    moshiDead = justAlleles[10:16]   #these specify which pools we are 
    tarime = justAlleles[16:22]
    #merus = justAlleles[22]
    Astack=[]
    Mstack=[]
    Tstack=[]
    
    for indv in moshiAlive:
        SNP = indv.split(':')
        numbers_int = [int(x) for x in SNP]
        total = sum(numbers_int)
        Astack.append(total)
    
    for indv in moshiDead:
        SNP = indv.split(':')
        numbers_int = [int(x) for x in SNP]
        total = sum(numbers_int)
        Mstack.append(total)
        
    for indv in tarime:
        SNP = indv.split(':')
        numbers_int = [int(x) for x in SNP]
        total = sum(numbers_int)
        Tstack.append(total)  
  
    for nums in Astack:
        stack_int = [int(x) for x in Astack]    
        Astack_tot = sum(stack_int) / len(stack_int)
    
    for nums in Mstack:
        stack_int = [int(x) for x in Mstack]    
        Mstack_tot = sum(stack_int) / len(stack_int)
        
    for nums in Tstack:
        stack_int = [int(x) for x in Tstack]    
        Tstack_tot = sum(stack_int) / len(stack_int)
        
    #print (chrom,"\t",where,"\t",Astack_tot,"\t",Mstack_tot,"\t",Tstack_tot)
    out.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, where,Astack_tot, Mstack_tot, Tstack_tot))
    
sync.close()
out.close()

depth = open("MoshiA_MoshiD_Tarime_depth_2R_HEAD", "r")
out = open('MoshiD_Tarime_winDepth_2_2R_HEAD', 'w')
#N=2
#for i in range(N): #check head of file to make sure correct
#    line=depth.readline().strip()
#    print(line)

counter = 0
winsize = 2 #set window depth here
Mstack=[]
Tstack=[]
Astack=[]

#print("chr\tstart\tend\tmid\tMoshiDead\tTarime")
out.write("chr\tstart\tend\tmid\tMoshiAlive\tMoshiDead\tTarime\n")

for position in depth:
    position = position.rstrip()
    if "chr" not in position:
        counter +=1
        
        if counter < winsize:
            lineSplit = position.split()
            if counter == 1:
                start = lineSplit[1]
            Astack.append(lineSplit[2])
            Mstack.append(lineSplit[3])
            Tstack.append(lineSplit[4])

        if counter == winsize:
            lineSplit = position.split()
            contig = lineSplit[0]
            end = lineSplit[1]
            Astack.append(lineSplit[2])
            Mstack.append(lineSplit[3])
            Tstack.append(lineSplit[4])
            for nums in Astack:
                Astack_int = [float(x) for x in Astack]    
                Astack_tot = sum(Astack_int) / len(Astack_int)
            for nums in Mstack:
                Mstack_int = [float(x) for x in Mstack]    
                Mstack_tot = sum(Mstack_int) / len(Mstack_int)           
            for nums in Tstack:
                Tstack_int = [float(x) for x in Tstack]    
                Tstack_tot = sum(Tstack_int) / len(Tstack_int)
                
            start = int(start)
            end = int(end)
            plus = (start + end) / 2
            #print(contig,"\t",start,"\t",end,"\t",plus,"\t",Mstack_tot,"\t",Tstack_tot)
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contig, start, end, plus, Astack_tot, Mstack_tot, Tstack_tot))
            del Astack[:]
            del Mstack[:]
            del Tstack[:]
            counter = 0

depth.close()
out.close()


import matplotlib.pyplot as plt
%matplotlib inline

out = open('MoshiA_MoshiD_Tarime_winDepth_1000_2R', 'r')

Position = []
Alive = []
Dead = []
Tarime = []


for line in out:
    line = line.rstrip()
    if "chr" not in line:
        lineSplit = line.split()
        Position.append(lineSplit[3])
        Alive.append(lineSplit[4])
        Dead.append(lineSplit[5])
        Tarime.append(lineSplit[6])
        
        
plt.scatter(Position, Alive)
