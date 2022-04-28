*************************************************************************
*     This program is to identify the best alignment of two protein 
*     structures to give the best TM-score. By default, TM-score is 
*     normalized by the second protein. The program can be freely 
*     copied or modified or redistributed.
*
*     On 2005/06/01, A small bug of two-point superposition was fixed.
*
*     On 2005/10/19, the program was reformed so that the alignment
*     results are not dependent on the specific compilers.
*
*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
*     (For comments, please email to: yzhang@ku.edu)
*************************************************************************

*************************************************************************
*
* Attribution Notice: 
*     This file comes from Yang Zhang's <yzhang@ku.edu> TMalign.f.
*     Yang Zhang's TMalign.f is copyrighted by Yang Zang
*     and he has all rights for that.
*     2006.07.04 Sunjoong LEE <sunjoong@gmail.com>
*         "Program compares" has been separated from subroutins
*         and modified to "subroutine TMalign."
*
*************************************************************************
      
      SUBROUTINE TMalign

      PARAMETER(nmax=3000)
      PARAMETER(nmax2=6000)
      PARAMETER(nmax3=9000)

      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002

      character*3 ss1(nmax),ss2(nmax)
      character seq1(0:nmax),seq2(0:nmax)
      character aseq1(nmax2),aseq2(nmax2),aseq3(nmax2)

      dimension m1(nmax),m2(nmax)
      dimension mm1(nmax),mm2(nmax)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/init/invmap_i(nmax)

      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d8/d8

      real rmsd8_al, TM8, seq_id
      integer n8_al

      common/result/rmsd8_al, TM8, seq_id, aseq1, aseq2, aseq3, n8_al

      common/pdb/ss1, ss2, mm1, mm2, seq1, seq2

      common/filenames/pdb1, pdb2, outname
      character*100 pdb1, pdb2, outname

      common/options/m_out, m_fix, m_ave, m_d0_min, m_d0,
     &               L_fix, d0_min_input, d0_fix
      integer m_out, m_fix, m_ave, m_d0_min, m_d0
      data m_out, m_fix, m_ave, m_d0_min, m_d0 / -1, -1, -1, -1, -1 /

*!!!  Scale of TM-score in search is based on the smaller protein --------->
      d0_min=0.5
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo
      
***** do alignment **************************
      CALL super_align          !to find invmap(j)
      
************************************************************
***   resuperpose to find residues of dis<d8 ------------------------>
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
      call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !TM-score with dis<d8 only

*!!!  Output TM-score is based on the second protein------------------>
      d0_min=0.5                !for output
      anseq=nseq2               !length for defining final TMscore
      if(m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 !<L>
      if(anseq.lt.anseq_min)anseq=anseq_min
      if(m_fix.eq.1)anseq=L_fix !input length
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      
***   remove dis>d8 in normal TM-score calculation for final report----->
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(ss1(m1(i)).eq.ss2(m2(i)))then
               n_eq=n_eq+1
            endif
         endif
      enddo
      seq_id=float(n_eq)/(n_al+0.00000001)
      n8_al=j
      d0_input=d0
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8*n8_al/anseq       !TM-score after cutoff
      
********* for output superposition ******************************
      if(m_out.eq.1)then
 1237    format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238    format('TER')
 1239    format('CONECT',I5,I5)
 900     format(A)
 901     format('select atomno=',I4)
         OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln
***   script:
         write(7,900)'load inline'
         write(7,900)'select atomno<2000'
         write(7,900)'wireframe .45'
         write(7,900)'select none'
         write(7,900)'select atomno>2000'
         write(7,900)'wireframe .20'
         write(7,900)'color white'
         do i=1,n8_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &           (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
            if(dis2.le.5)then
               write(7,901)m1(i)
               write(7,900)'color red'
               write(7,901)2000+m2(i)
               write(7,900)'color red'
            endif
         enddo
         write(7,900)'select all'
         write(7,900)'exit'
         write(7,104)pdb1,nseq1
 104     format('REMARK Chain 1:',A10,'  Size=',I4)
         write(7,105)pdb2,nseq2,int(anseq)
 105     format('REMARK Chain 2:',A10,'  Size=',I4,
     &        ' (TM-score is normalized by ',I4,')')
         write(7,106)n8_al,rmsd8_al,TM8,seq_id
 106     format('REMARK Aligned length=',I4,', RMSD=',f6.2,
     &        ', TM-score=',f7.5,', ID=',f5.3)
***   chain1:
         do i=1,n8_al
            write(7,1237)m1(i),ss1(m1(i)),mm1(m1(i)),
     &           xtm1(i),ytm1(i),ztm1(i)
         enddo
         write(7,1238)          !TER
         do i=2,n8_al
            write(7,1239)m1(i-1),m1(i) !connect atoms
         enddo
***   chain2:
         do i=1,n8_al
            write(7,1237)2000+m2(i),ss2(m2(i)),mm2(m2(i)),
     $           xtm2(i),ytm2(i),ztm2(i)
         enddo
         write(7,1238)
         do i=2,n8_al
            write(7,1239)2000+m2(i-1),2000+m2(i)
         enddo
      endif
*^^^^^^^^^^^^^^^^^^ output finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

************  output aligned sequences **************************
      ii=0
      i1_old=1
      i2_old=1
      do i=1,n8_al
         do j=i1_old,m1(i)-1
            ii=ii+1
            aseq1(ii)=seq1(j)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do j=i2_old,m2(i)-1
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seq2(j)
            aseq3(ii)=' '
         enddo
         ii=ii+1
         aseq1(ii)=seq1(m1(i))
         aseq2(ii)=seq2(m2(i))
         dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &     (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
         if(dis2.le.5)then
           aseq3(ii)=':'
         else
           aseq3(ii)='.'
         endif
         i1_old=m1(i)+1
         i2_old=m2(i)+1
      enddo
      do i=i1_old,nseq1
         ii=ii+1
         aseq1(ii)=seq1(i)
         aseq2(ii)='-'
         aseq3(ii)=' '
      enddo
      do i=i2_old,nseq2
         ii=ii+1
         aseq1(ii)='-'
         aseq2(ii)=seq2(i)
         aseq3(ii)=' '
      enddo

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

