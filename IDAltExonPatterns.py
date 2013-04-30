"""
Joshua Arribere April 20, 2013 (woooo 4-20)

Input: ucscAnnots.txt - ucsc-formatted annotations, i.e.
        #name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds proteinID alignID
        except that the spaces are tabs

Output: The number of different alternative splicing events, classified into
    different events

run as python IDAltExons.py ucscAnnots.txt outPrefix
April 21, 2013 - JOSH decided to go a different route. Rather than make a txt
    graph for each txt, I decided to simply compare the pattern of splice sites.
April 23, 2013 - JOSH changed to do a slightly different way--decision tree.
April 28, 2013 - JOSH changes how AFEs are called. Just requires two FEs that
    are non-overlapping
April 29, 2013 - JOSH copied over from IDAltExonPatterns5.py and uploaded to git
"""
import sys, common, csv, numpy, collections
from pyx import *

def parseAnnots(annots):
    """Will parse annotations to a dict of gene:info"""
    a={}
    with open(annots,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        freader.next()#skips the header
        for row in freader:
            gene=row[0]
            chr=row[1]
            strand=row[2]
            exonStarts=map(int, row[8].strip(',').split(','))
            exonEnds=map(int, row[9].strip(',').split(','))
            cdsStart=int(row[5])
            cdsEnd=int(row[6])
            txtStart=int(row[3])
            txtEnd=int(row[4])
            
            TL=1
            if strand=='+':
                if txtStart==cdsStart:
                    TL=0
            else:
                if txtEnd==cdsEnd:
                    TL=0
            
            if cdsStart!=cdsEnd and TL==1:#then it's a coding gene
                a[gene]={'chr':chr,'strand':strand,
                         'exonStarts':exonStarts,'exonEnds':exonEnds,
                         'cdsStart':[cdsStart],'cdsEnd':[cdsEnd],
                         'txtStart':[txtStart],'txtEnd':[txtEnd]}
    
    return a

def getGroups(b,txt,done_txts):
    """Will go through all keys in b[txt] and check that their keys are already
    in done_txts. If not, will call this function with the new txt"""
    if txt not in done_txts:
        done_txts.append(txt)
    for txt2 in b[txt]:
        if txt2 not in done_txts:
            done_txts=getGroups(b,txt2,done_txts)
        else:
            pass
    return done_txts

def getTxtsThatShareSS(annots):
    """After talking to Merkin, he says that for grouping transcripts together
    it's better to group by those that share at least one splice site.
    Previously I was just going by txts that overlapped at all. And Jason said
    there are a fair number of txts that overlap UTRs (e.g. 5'UTR overlaps w/
    3'end of another gene)."""
    #First step is to make a map of the genome with all coding exon positions
    #telling which txts they come from.
    a={'+':collections.defaultdict(lambda:collections.defaultdict(list)),
       '-':collections.defaultdict(lambda:collections.defaultdict(list))}
    
    for txt in annots:
        chr=annots[txt]['chr']
        strand=annots[txt]['strand']
        
        exons=zip(annots[txt]['exonStarts'],annots[txt]['exonEnds'])
        
        for exon in exons:
            for ss in exon:
                a[strand][chr][ss].append(txt)
    
    #Now will make a txt dictionary where the keys are txt1 IDs, values are
    #dictionaries of txt2:1 where a 1 indicates that txt shares a coding region
    #with txt1
    b=collections.defaultdict(dict)
    for strand in a:
        for chr in a[strand]:
            for position in a[strand][chr]:
                txts=a[strand][chr][position]
                
                for txt1 in txts:
                    for txt2 in txts:
                        b[txt1][txt2]=1
    
    c=[]
    done={}
    for txt in b:
        if txt not in done:
            group=getGroups(b,txt,[txt])
            c.append(group)
            
            for txt2 in group:
                done[txt2]=1
        else:
            pass
    
    return c

def mkPattern2(group,annots):
    """Given group, a list of txts, will return the pattern of start/ss5/ss3/end
    that makes up the txt as a dictionary of {txt:{position:start/ss5/ss3/end}}
    """
    strand=annots[group[0]]['strand']
    
    a=[]
    if strand=='+':
        positions=list(set([entry for txt in group for key in ('exonStarts','exonEnds',
                                                      'txtStart')
                   for entry in annots[txt][key]]))
        positions.sort()
        for position in positions:
            a.append([])
            for txt in group:
                if position in annots[txt]['txtStart']:
                    a[-1].append('txtStart')
                elif position in annots[txt]['exonStarts']:
                    a[-1].append('ss3')
                elif position in annots[txt]['exonEnds']:
                    a[-1].append('ss5')
                #elif position in annots[txt]['cdsStart']:
                #    a[-1].append('cdsStart')
                else:
                    a[-1].append(0)
    
    elif strand=='-':
        positions=list(set([entry for txt in group for key in ('exonStarts','exonEnds',
                                                      'txtEnd')
                   for entry in annots[txt][key]]))
        positions.sort()
        positions.reverse()
        for position in positions:
            a.append([])
            for txt in group:
                if position in annots[txt]['txtEnd']:
                    a[-1].append('txtStart')
                elif position in annots[txt]['exonStarts']:
                    a[-1].append('ss5')
                elif position in annots[txt]['exonEnds']:
                    a[-1].append('ss3')
                #elif position in annots[txt]['cdsEnd']:
                #    a[-1].append('cdsStart')
                else:
                    a[-1].append(0)
    
    else:
        print 'Strand not recognized.', sys.exit()
    
    return positions,a

def combinePatterns(pattern1,pattern2,strand):
    """patterni is a dict of mkPattern[txt]. Will combine the patterns into a
    list of txtStart,ss5,ss3,txtEnd,cdsStart."""
    positions=list(set(pattern1.keys()).union(pattern2.keys()))
    positions.sort()
    if strand=='-':
        positions.reverse()
    pattern=[]
    for position in positions:
        if (position in pattern1) and (position in pattern2):
            pattern.append([pattern1[position],pattern2[position]])
        elif position in pattern1:
            pattern.append([pattern1[position]])
        elif position in pattern2:
            pattern.append([pattern2[position]])
        else:
            print 'A position was randomly generated.', sys.exit()
    
    return pattern

class txtPattern(list):
    """Simple class that will be used to easily loop through a pattern list and
    find the next non-zero entry of a list"""
    
    def iterate(pattern):
        """Generator that will loop through all the values of pattern and return
        those that are not (0,0). Will return (index,pattern[index])"""
        for ii in range(len(pattern)):
            if pattern[ii]!=(0,0):
                yield ii,pattern[ii]
    
    def next(self,jj):
        """Will find the next non-(0,0) entry"""
        jj+=1
        if jj>=len(self):
            return (0,0)#gone beyond boundary
        if self[jj]!=(0,0):
            return jj,self[jj]
        else:
            return self.next(jj)
    
    def prev(self,jj,kk):
        """Will start with the kkth entry of the jjth entry of self (self[jj][kk]
        and return the first preceding non-zero instance, walking up one entry
        at a time in self. If jj becomes negative, will return 0"""
        jj-=1
        if jj==-1:
            return 0
        if self[jj][kk]!=0:
            return self[jj][kk]
        else:
            return self.prev(jj,kk)
    

def getcdsStartIndex(group,positions,annots):
    """Will go through each txt in group, ID the cdsStart, and return a list of
    indexes of the first upstream position in positions"""
    strand=annots[group[0]]['strand']
    a=[]
    if strand=='+':
        for txt in group:
            cdsStart=annots[txt]['cdsStart'][0]
            for i in range(len(positions)):
                if cdsStart-positions[i]<0:
                    a.append(i-1)#go one previous since this one is downstream of cdsStart (negative)
                    break
    elif strand=='-':
        for txt in group:
            cdsStart=annots[txt]['cdsEnd'][0]
            for i in range(len(positions)):
                if cdsStart-positions[i]>=0:
                    a.append(i)
                    break
    return a

def extendUpstream(FEIndexes,txtStart):
    """Given FEIndexes, which is a list of lists of [txtStart,ss5] positions and
    txtStart, will find the leftmost boundary of all overlapping tuples"""
    for tuple in FEIndexes:
        if tuple[1]>=txtStart:
            if tuple[0]<txtStart:
                return extendUpstream(FEIndexes,tuple[0])
    return txtStart

def extendDownstream(FEIndexes,ss5):
    """Given FEIndexes, which is a list of lists of [txtStart,ss5] positions and
    ss5, will find the rightmost boundary of all overlapping tuples"""
    for tuple in FEIndexes:
        if tuple[0]<ss5:
            if tuple[1]>ss5:
                return extendDownstream(FEIndexes,tuple[1])
    return ss5

def getFirstExons(group,patterns):
    """Will return a list of the First Exons of the txts in groups as index in
    patterns. EDIT: will attempt to extend exons to encompass all overlapping
    FEs into a single FE"""
    FEIndexes=[]#will keep track of the indexes of the first exon, i.e.
    #[[0,4],[5,6],[5,7]] might be the pattern for a gene with an alt 5'SS and
    #2 AFEs
    for ii in range(len(group)):
        for jj in range(len(patterns)):
            if patterns[jj][ii]=='txtStart':
                FEIndexes.append([jj])
            elif patterns[jj][ii]=='ss5':
                FEIndexes[-1].append(jj)
                break
    
    assert len(FEIndexes)==len(group), 'A txt was lost or generated.'
    FEIndexesOverlap=[]
    for ii in range(len(FEIndexes)):
        txtStartOverlap=extendUpstream(FEIndexes,FEIndexes[ii][0])
        ss5Overlap=extendDownstream(FEIndexes,FEIndexes[ii][1])
        if [txtStartOverlap,ss5Overlap] not in FEIndexesOverlap:
            FEIndexesOverlap.append((txtStartOverlap,ss5Overlap))
        
    return FEIndexesOverlap

def nonOverlappingAndUpstream(FEs_i,FEs_j):
    """Will figure out with these two tuples are non-overlapping and FEs_i is
    upstream of FEs_j"""
    if FEs_j[0]>FEs_i[1]:
        return 1
    else:
        return 0

def getFEsWithStarts(cdsStarts,patterns):
    """Will return a list of all non-redundant first exons with a start codon"""
    a=[]
    for ii in range(len(cdsStarts)):
        for jj in range(len(patterns)):
            tuple=[]
            if patterns[jj][ii]=='txtStart':
                tuple.append(jj)
                for kk in range(jj+1,len(patterns)):
                    if patterns[kk][ii]=='ss5':
                        if cdsStarts[ii]<=kk:
                            if tuple not in a:
                                a.append(tuple+[kk])
                                break
                        else:
                            break
                break
    return a

def countClasses(groups,annots):
    """groups is a list of lists with txt IDs, whose annotation information is
    given in annots[txt]"""
    
    metaDone=collections.defaultdict(int)
    metaAllDone=collections.defaultdict(int)
    for group in groups:
        positions,patterns=mkPattern2(group,annots)
        done=collections.defaultdict(dict)#will need to keep track so don't count the same event twice
        allDone=collections.defaultdict(dict)#will be used to keep track of all alternative events, whether TL affecting or not
        cdsStarts=getcdsStartIndex(group,positions,annots)
        
        FEs=getFirstExons(group,patterns)
        allDone['AFE']=dict((str(k),1) for k in FEs[:-1])
        done['AFE']=dict((str(k),1) for k in FEs[:-1])
        #This does not count the number of AFEs, but rather counts the pairs. One pair is one event, not two exons. This is done to be similar to the other alternative events.
        
        FEStarts=getFEsWithStarts(cdsStarts,patterns)
        done['AFEStart']=dict((str(k),1) for k in FEStarts)
        
        for i,j in [(i,j) for i in range(len(group)) for j in range(len(group)) if i!=j]:
            txti,txtj=group[i],group[j]
            patternij=txtPattern([(entry[i],entry[j]) for entry in patterns])
            #cdsStarts=[patternij.index(entry) for entry in patternij if 'cdsStart' in entry]
            cdsStartsij=(cdsStarts[i],cdsStarts[j])
            
            for ii,event in patternij.iterate():#thus begins the decision tree
                if event==('ss5','ss5'):
                    jj,nextevent=patternij.next(ii)
                    if nextevent==('ss3',0):
                        altExon1=jj#defines the first skipped exon as the alternative one
                        jj,nextevent=patternij.next(jj)
                        if nextevent==('ss5',0):
                            jj,nextevent=patternij.next(jj)
                            if nextevent==(0,'ss3'):
                                altExon2=jj
                                jj,nextevent=patternij.next(jj)
                                if nextevent==(0,'ss5'):
                                    jj,nextevent=patternij.next(jj)
                                    if nextevent==('ss3','ss3'):
                                        #Now know it's a Mutually Exclusive Codon. Let's ID the cdsStart locations and figure out exactly what type
                                        allDone['ME'][str(altExon1)+str(altExon2)]=1
                                        if min(cdsStartsij)>=jj:#then the start is downstream of all this alternative stuff
                                            done['ME'][str(altExon1)+str(altExon2)]=1#example: uc007tih.2 uc007tii.2
                                        elif cdsStartsij[0]==altExon1 and cdsStartsij[1]==altExon2:
                                            done['MEStart2'][str(altExon1)+str(altExon2)]=1#positive example: uc008xzl.1 uc012dxz.1
                                        elif max(cdsStartsij)>=jj:
                                            if cdsStartsij[0]==altExon1 or cdsStartsij[1]==altExon2:
                                                done['MEStart1'][str(altExon1)+str(altExon2)]=1#positive example: uc008riy.1 uc012cxu.1
                            elif nextevent==('ss3','ss3'):
                                #It's got the architecture of a Skipped Exon. ID the location of cdsStarts next
                                allDone['SE'][altExon1]=1
                                if min(cdsStartsij)>=jj:
                                    done['SE'][altExon1]=1#example: uc009elk.1 uc009ell.1
                                elif cdsStartsij[0]==altExon1 and cdsStartsij[1]>=jj:
                                    done['SEStart'][altExon1]=1#example: uc008brw.1 uc008brx.1
                        elif nextevent==(0,'ss3'):
                            altExon2=jj
                            #Then it's a alternative 3'ss. Let's ID the cdsStart locations
                            allDone['SS3'][str(altExon1)+str(altExon2)]=1
                            if min(cdsStartsij)>=jj:#example: uc007bso.2 uc011wob.1
                                done['SS3'][str(altExon1)+str(altExon2)]=1#keep track of the alternative event as a pair. That way if I detect another alternative event with different bounds I'll count it again.
                            elif cdsStartsij[0]==altExon1:#xample uc007nfa.1 uc007nfb.2
                                done['SS3Start'][str(altExon1)+str(altExon2)]=1
                elif event==('ss5',0):
                    jj,nextevent=patternij.next(ii)
                    if nextevent==(0,'ss5'):
                        altExon1=jj
                        jj,nextevent=patternij.next(jj)
                        if nextevent==('ss3','ss3'):
                            altExon2=jj
                            #Then it's an alternative 5'ss. Where are the cdsStarts?
                            allDone['SS5'][str(altExon1)+str(altExon2)]=1
                            if min(cdsStartsij)>=jj:#example uc007hkz.1 uc007hky.1
                                done['SS5'][str(altExon1)+str(altExon2)]=1#another example: uc008sqn.1 uc008sqo.1
                            elif cdsStartsij[1]==altExon1:
                                done['SS5Start'][str(altExon1)+str(altExon2)]=1#example uc007hba.1 uc007hbb.1
                    elif nextevent==(0,'txtStart'):
                        altExon2=jj
                        jj,nextevent=patternij.next(jj)
                        #if nextevent==(0,'ss5'):
                        #    jj,nextevent=patternij.next(jj)
                        #    if nextevent==('ss3','ss3'):
                        #        #Then it's an Alternative First Exon Pattern. Let's ID where the cdsStarts lie in relation.
                        #        allDone['AFE'][str(ii)+str(altExon2)]=1
                        #        if min(cdsStartsij)>=jj:
                        #            done['AFE'][str(ii)+str(altExon2)]=1#example uc007qly.1 uc007qlz.1
                        #        elif max(cdsStartsij)>=jj:
                        #            if cdsStartsij[0]<ii or cdsStartsij[1]==altExon2:#note this requires the cdsStart of the longer isoform to be upstream of the alternative event
                        #                done['AFEStart1'][str(ii)+str(altExon2)]=1#example uc012bsj.1 uc008iuq.2
                        #        elif cdsStartsij[0]<ii and cdsStartsij[1]==altExon2:
                        #            done['AFEStart2'][str(ii)+str(altExon2)]=1#example uc007qma.1 uc007qmb.1
                        if nextevent==('ss3',0):
                            #then it looks like an Alternative First Exon spliced internal to another exon that's also the Txn start site for the other isoform. AFEASS. Find cdsStarts
                            allDone['SFiEASS'][str(altExon2)+str(jj)]=1
                            if min(cdsStartsij)>=jj:#example uc009aek.1 uc009ael.1
                                done['SFiEASS'][str(altExon2)+str(jj)]=1
                            elif max(cdsStartsij)>=jj:
                                if cdsStartsij[1]==altExon2:#example uc008qgm.1 uc008qgp.1
                                    done['SFiEASSStart1'][str(altExon2)+str(jj)]=1#cdsStart between txtStart and 3'ss
                                elif cdsStartsij[0]<ii:#two positive pair confirmed: uc009kpa.1 uc009koy.1; uc009kpa.1 uc009koz.1
                                    done['SFiEASSStart2'][str(altExon2)+str(jj)]=1#cdsStart upstream of txtStart
                            elif cdsStartsij[1]==altExon2 and cdsStartsij[0]<ii:#example uc009jzo.2 uc009jzp.1
                                done['SFiEASSStart3'][str(altExon2)+str(jj)]=1#cdsStart upstream of txtStart; second cdsStart between txtStart and 3'ss
                    elif nextevent==('ss3','txtStart'):
                        #Then an isoform is spliced to the txtStart site of another isoform. Figure out where the cdsStarts are.
                        allDone['SFE'][jj]=1
                        if min(cdsStartsij)>=jj:#example uc008nrn.1 uc008nrt.1
                            done['SFE'][jj]=1
                        elif max(cdsStartsij)>=jj:
                            if cdsStartsij[0]<ii:#positive example: uc007rlf.1 uc007rlg.1
                                done['SFEStart'][jj]=1
                    else:
                        prevevent=patternij.prev(ii,1)
                        if prevevent in ('ss3','txtStart'):
                            jj,nextevent=patternij.next(ii)
                            if nextevent==('ss3',0):
                                #great, it looks like a retained intron. Where are the cdsStarts in relation?
                                allDone['IR'][str(ii)+str(jj)]=1
                                if min(cdsStartsij)>=jj:
                                    done['IR'][str(ii)+str(jj)]=1#example uc007gae.1 uc007gaf.1
                                elif max(cdsStartsij)>=jj and cdsStartsij[1]==ii:
                                    done['IRStart'][str(ii)+str(jj)]=1#example uc009gfd.1 uc009gfb.1
                elif event==('txtStart',0):
                    jj,nextevent=patternij.next(ii)
                    if nextevent==(0,'txtStart'):
                        #We're calling that alternative txtStart sites--they are consecutive with no other events between them. Where do the cdsStarts lie in relation?
                        allDone['ASS'][str(ii)+str(jj)]=1
                        if min(cdsStartsij)>=jj:#positive example: uc007hkz.1 uc007hky.1
                            done['ASS'][str(ii)+str(jj)]=1
                        elif cdsStartsij[0]==ii:#example uc009feu.2 uc009fex.2
                            done['ASSStart'][str(ii)+str(jj)]=1
        for key in done:
            metaDone[key]+=len(done[key])
        for key in allDone:
            metaAllDone[key]+=len(allDone[key])
    
    
    for key in metaDone:
        if key in metaAllDone:
            #print key,'\t%.4f'%(metaDone[key]/float(sum(metaDone.values()))),'\t%.4f'%(metaAllDone[key]/float(sum(metaAllDone.values())))
            print key+'\t'+str(metaDone[key])+'\t'+str(metaAllDone[key])
        else:
            #print key,'\t%.4f'%(metaDone[key]/float(sum(metaDone.values())))
            print key+'\t'+str(metaDone[key])
    print sum(metaDone.values())
    print sum(metaAllDone.values())
    sys.exit()
    return metaDone

def mkTable(classCounts,outPrefix):
    """This will make a table summarizing the abundances of different events"""
    print classCounts.keys()
    print len(classCounts)
    c=canvas.canvas()
    text.set(mode="latex")
    types=['SE',
           'SEStart',
           'SS5',
           'SS5Start',
           'SS3',
           'SS3Start',
           'AFE',
           'AFEStart1',
           'AFEStart2',
           'ME',
           'MEStart1',
           'MEStart2',
           'ASS',
           'ASSStart',
           'SFE',
           'SFEStart',
           'SFiEASS',
           'SFiEASSStart1',
           'SFiEASSStart2',
           'SFiEASSStart3',
           'IR',
           'IRStart']
    #There are 9 broad classes, and 2-4 members within each class. I think this would be
    #well represented as a table with members aligned vertically
    broadTypes=['SE',
           'SS5',
           'SS3',
           'AFE',
           'ME',
           'ASS',
           'SFE',
           'SFiEASS',
           'IR']
    totals=dict((key,0) for key in broadTypes)
    for key1 in broadTypes:
        for key2 in classCounts:
            if key2.startswith(key1):
                totals[key1]+=classCounts[key2]
    #first set parameters that will be useful for easy, general changes to the
    #table
    hSpacing=10.
    vSpacing=4.
    textWidth=8.
    #defaultAlign=[text.parbox(textWidth),text.halign.boxcenter,text.valign.bottom]
    defaultAlign=[text.parbox(textWidth, baseline=text.parbox.bottom),
                  text.halign.boxcenter,text.halign.flushcenter]
    titleAlign=[text.parbox(textWidth),text.valign.middle,
                  text.halign.boxright,text.halign.flushright,
                  text.size(2)]
    titleAlign2=[text.parbox(textWidth, baseline=text.parbox.bottom),
                  text.halign.boxcenter,text.halign.flushcenter,
                  text.size(4)]
    breakdownAlign=[text.parbox(10.),text.valign.middle,
                  text.halign.boxleft,text.halign.flushleft,
                  text.size(2)]
    diagramWidth=6.
    exonWidth=diagramWidth/4.
    exonHeight=vSpacing/4.
    lineWidth=[style.linewidth.normal]
    
    
    #line across the top
    c.stroke(path.line(-textWidth/2.,-vSpacing/8.,diagramWidth*2+hSpacing,-vSpacing/8.),
             [style.linewidth.THICK])
    #write column headers across the top
    c.text(0,0,'Event Type (Number)',titleAlign2)
    c.text(diagramWidth/2.+hSpacing/2.,0,'Schematic',titleAlign2)
    c.text(diagramWidth*2+hSpacing/2.,0,'Start Site Location Breakdown',titleAlign2)
    
    
    
    #Skipped Exon
    ii=0
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Skipped Exon ('+str(totals['SE'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['SE'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['SEStart']),breakdownAlign)
    #draw a left exon
    c.fill(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #and the alternative exon SJs
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*1.25,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1.25,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*2.75,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.75,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    
    
    #SS5
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Alternative 5\'SS ('+str(totals['SS5'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['SS5'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['SS5Start']),breakdownAlign)
    #draw a left exon
    c.fill(path.rect(-hSpacing/2.+hSpacing,-ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             [style.linewidth.normal])
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*1,-ii*vSpacing-5*vSpacing/8., exonWidth/2., exonHeight),
           [color.rgb.black])
    #and the alternative exon SJs
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*(1.5+(3-1.5)/2.),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(1.5+(3-1.5)/2.),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    
    #SS3
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Alternative 3\'SS ('+str(totals['SS3'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['SS3'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['SS3Start']),breakdownAlign)
    #draw a left exon
    c.fill(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             [style.linewidth.normal])
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*2.5, -ii*vSpacing-5*vSpacing/8., exonWidth/2., exonHeight),
           [color.rgb.black])
    #and the alternative exon SJs
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*(1+(2.5-1)/2.),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(1+(2.5-1)/2.),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    
    #AFE
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Alternative First Exon ('+str(totals['AFE'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['AFE'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['AFEStart']),breakdownAlign)
    #draw a left exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*1.5, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #draw the txtStarts
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*0,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*0,-ii*vSpacing-4*vSpacing/8.+exonHeight),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*0,-ii*vSpacing-4*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*0.25,-ii*vSpacing-4*vSpacing/8.+exonHeight),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-6*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1.5,-ii*vSpacing-6*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*1.75,-ii*vSpacing-6*vSpacing/8.),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    #and the alternative exon SJs
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*2.75,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.75,-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    
    #ME
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Mutually Exclusive Exons ('+str(totals['ME'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['ME'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['MEStart1'])+'\n\nBoth Start Sites in Alt Exons: '+
           str(classCounts['MEStart2']),breakdownAlign)
    #draw a left exon
    c.fill(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw the alternative exons
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*(4/3.), -ii*vSpacing-5*vSpacing/8., exonWidth/2., exonHeight),
           lineWidth)
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*(3-5./6), -ii*vSpacing-5*vSpacing/8., exonWidth/2., exonHeight),
           lineWidth)
    #draw splicing of the top isoform
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*(7./6),-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(7./6),-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*(4./3),-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(1./2+4./3),-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*(3-5./6+1./4),-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(3-5./6+1./4),-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the splicing of the bottom isoform
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*(1+(2-5./6)/2),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(1+(2-5./6)/2),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*(1+(2-5./6)),-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(3-2./6),-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*(3-1./6),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*(3-1./6),-ii*vSpacing-5*vSpacing/8.-exonHeight*2./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.),
             lineWidth)
    
    #ASS
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Alternative Transcription Start Site ('+str(totals['ASS'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['ASS'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['ASSStart']),breakdownAlign)
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*2, -ii*vSpacing-5*vSpacing/8., exonWidth*2, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*1, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #draw the txtStarts
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*1,-ii*vSpacing-6*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*1,-ii*vSpacing-6*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*1.25,-ii*vSpacing-6*vSpacing/8.),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-4*vSpacing/8.+exonHeight),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-4*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2.25,-ii*vSpacing-4*vSpacing/8.+exonHeight),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    
    #SFE
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Skipped First Exon ('+str(totals['SFE'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['SFE'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['SFEStart']),breakdownAlign)
    #draw the left exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the txtStart
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-6*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-6*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-6*vSpacing/8.),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    
    #SFiEASS
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Skipped First Exon and Alternative Transcription Start Site ('+str(totals['SFiEASS'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['SFiEASS'])+'\n\nStart Site in First Alt Exon: '+
           str(classCounts['SFiEASSStart2'])+'\n\nStart Site in Second Alt Exon: '+
           str(classCounts['SFiEASSStart1'])+'\n\nStart Site in Both Alt Exons: '+
           str(classCounts['SFiEASSStart3']),breakdownAlign)
    #draw the left exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           lineWidth)
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the txtStart
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-5*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-6*vSpacing/8.),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2.5,-ii*vSpacing-6*vSpacing/8.,
                       -hSpacing/2.+hSpacing+exonWidth*2.75,-ii*vSpacing-6*vSpacing/8.),
        [style.linewidth.normal, deco.earrow([deco.stroked([color.rgb.black]),
                       deco.filled([color.rgb.black])], size=exonWidth/8.)])
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*2.5, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    
    #IR
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing-vSpacing/2.,'Intron Retention ('+str(totals['IR'])+')',titleAlign)
    #breakdown by start sites and do counts
    c.text(diagramWidth+3*hSpacing/4.,-ii*vSpacing-vSpacing/2.,
           'Start Site Downstream: '+str(classCounts['IR'])+'\n\nStart Site in Alt Exon: '+
           str(classCounts['IRStart']),breakdownAlign)
    #draw a left exon
    c.fill(path.rect(-hSpacing/2.+hSpacing, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw a start over the right exon
    c.text(-hSpacing/2.+hSpacing+exonWidth*3.5,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.5,
           'Start',[color.cmyk(0.97,0,0.75,0)]+defaultAlign)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*3.25,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25,
                       -hSpacing/2.+hSpacing+exonWidth*3.75,-ii*vSpacing-5*vSpacing/8.+exonHeight*1.25),
        [color.cmyk(0.97,0,0.75,0),
         style.linewidth.normal, deco.earrow([deco.stroked([color.cmyk(0.97,0,0.75,0)]),
                       deco.filled([color.cmyk(0.97,0,0.75,0)])], size=exonWidth/8.)])
    #draw a right exon
    c.fill(path.rect(-hSpacing/2.+hSpacing+exonWidth*3, -ii*vSpacing-5*vSpacing/8., exonWidth, exonHeight),
           [color.rgb.black])
    #draw splicing between left and right exons
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth,-ii*vSpacing-5*vSpacing/8.+exonHeight,
                       -hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3),
             lineWidth)
    c.stroke(path.line(-hSpacing/2.+hSpacing+exonWidth*2,-ii*vSpacing-5*vSpacing/8.+exonHeight*5./3,
                       -hSpacing/2.+hSpacing+exonWidth*3,-ii*vSpacing-5*vSpacing/8.+exonHeight),
             lineWidth)
    #draw the alternative exon
    c.stroke(path.rect(-hSpacing/2.+hSpacing+exonWidth*1, -ii*vSpacing-5*vSpacing/8., exonWidth*2, exonHeight),
           lineWidth)
    
    
    #total on the bottom
    ii+=1
    c.text(hSpacing/3.,-ii*vSpacing,'Total Events: '+str(sum(totals.values())),
           [text.parbox(textWidth),text.valign.middle,
                  text.halign.boxright,text.halign.flushright,
                  text.size(4)])
    
    c.writePDFfile(outPrefix)

def main(args):
    annots,outPrefix=args[0:]
    
    annots=parseAnnots(annots)
    
    groups=getTxtsThatShareSS(annots)
    
    groups=[entry for entry in groups if len(entry)>1]#restricts to only those txt groups with >1 txt
    
    classCounts=countClasses(groups,annots)
    common.rePickle(classCounts,outPrefix+'.p')
    classCounts=common.unPickle(outPrefix+'.p')
    mkTable(classCounts,outPrefix)

if __name__=='__main__':
    main(sys.argv[1:])