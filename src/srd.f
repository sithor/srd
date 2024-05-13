      subroutine srd(d,f,n,q,x,y,e,c,crit,iseed,thin,doptm,dactm,skinxx)
c%%!dir$ ATTRIBUTES DLLEXPORT :: SRD
c%%!dir$ ATTRIBUTES C, REFERENCE, ALIAS : 'srd_' :: SRD
c
c    d=data matrix  n x q 
c    f=frequency vector length n (integer)
c    n= number of rows of data matrix
c    q=number of attributes (rectangles)
c    x  is output 4 x q  matrix of x rectangle coordinates      
c    y  is output 4 x q  matrix of y rectangle coordinates
c    e is % error
c    c is 5 x 2**q matrix of cells statistics
c         1,2  are x,y coordinates of midpoint
c         3  is cell frequency    
c         4  is cell area of fitted configuration
c         5  is cell error   

c   crit is criterion for fitting  1=least squares (default)
c    2 is least abs difference,  3 is  minus log-likelihood
c   iseed is seed for random number generation
c  thin  is thinness parameter (=skinpa)

       integer q, proxyq
      integer d(n,q)
	  real f(n)
      common/ngrp/proxyq

      double precision x(4,7), y(4,7) 
c 4=#corners of rectangles, 7 = <=6 rectangles + boundary

      common/grpf/tot(6)
      Real tot

      integer bin(6)
      real ndisj(0:1023)
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      common/frqn/no(0:63)
      common/area/xn(0:63)
      common/posn/xm(0:63),ym(0:63)
      Real no, xn, xnopt(0:63),xmopt(0:63),ymopt(0:63)
      common/option/skinpa,absent , skinx    
      common/object/lad,lsd,logl,chi2

      logical lad,lsd,logl,chi2, nudge
      common/space/whtesp,greysp,blcksp  
      double precision dopt, doptm, dactm, dact
      Real er,err
      character*9 inital, initx
      character*50 string
      integer iarray(2)
      Real xng(6), ng(6)

      Real e
      Real c(5,64)
      integer crit
      common/error/eperc
      common/actual/dact
c      iseed=123456
10000  continue  
      lsd=.false.
      lad=.false.
      chi2=.false.
      logl=.false. 
      xsamll=0.0
      ysmall=0.0
      dx=1.0
      dy=1.0
	  	    call intpr("crit",4,crit,1)
      if(crit.eq.1) then
         lsd=.true.
      end if
          
      if(crit.eq.2) then
         lad=.true.
      end if
 
      if(crit.eq.3) then
          logl=.true.
          end if
 
      if(crit.eq.4) then
          chi2=.true.
          end if

c      skinpa=0.05
c      absent=0.05
c       skinpa=thin
      absent=0.0
        
      
      proxyq=q  
c  proxyq =  a device to allow q to be used in func(p)
c  first collapse data into cell counts of 2^q  table      
      ndisj=0  
      note1=0
      do 145 i=1,n
          if(f(i).lt.0) note1=1
         do 146 j=1,q
          bin(j)=d(i,j)  
146       continue
c  establish cell # k correponding to bin()         
         call ibinno(k,q,bin)
c     write(*,*)'k=', k
         ndisj(k)=ndisj(k)+f(i)
145    continue
      t=0
      xl=0
      do 144 k=0,2**q-1
          no(k)=ndisj(k)
          t=t+no(k)
          call binno(k,q,bin)
c     write(*,*)k, (bin(j),j=1,q), ndisj(k)
 
144   continue
       skinpa=thin
c  input thin <0 is same as thin=NA in srd.R, automatic penalising      
      if(thin.lt.0) then
	  
c  WHAT IS THIS????  why 0a fraction of sample size?? 
c   thinness penalising is  Dopt <- Dopt + max(eccentricity)*skinpa
          skinpa=0.05*t
      endif
      
      
      whtesp=ndisj(0)
      blcksp=ndisj(2**q-1)
      greysp=t-whtesp-blcksp
      
c      write(*,*)'spaces', whtesp, greysp,blcksp
c      set up inital configuration 

c      if(note1.eq.1) write(*,*)'Note:negative frequency count(s)'
c     if(logl) write(*,*) 'Fitting criterion D is -log-likelihood'
c      if(lsd)  write(*,*) 'Fitting criterion D is least squares'
c      if(lad)  write(*,*) 'Fitting criterion D is least abs. difference'
c      if(chi2)  write(*,*) 'Fitting criterion D is chi-square'

      call marg(q)
      doptm=1.0e20
c--------------------------------------------------------
c  loop over initial configurations 
c   input x,y. Set initial values from them. Attempt to formulate
c		"nudge"
		
		iq1 = q +1
		 do 876 i=1, iq1
		 lx(i) = x(1,i)
		 ly(i) = y(1,i)
		 rx(i) = x(3,i)
		 ry(i) = y(3,i)
876  	continue 
		nudge = .true. 
	
		if( isnan(lx(1)) ) nudge = .false. 
c		call realpr("lx(1)",5,lx(1),1)
c		call realpr("x(1,i)",6,x(1,i),1)
c		inital ="random"
c   ntries indicates inital configuration method >5 is random
      ntries=1
	  ntriex=5
	  if(nudge) ntriex =1
100   continue      

      if(.not.nudge) call gstart(q, ntries, inital,iseed)
      call wobble(q, Dopt)
cc         write(*,762)inital, ntries, Dopt,eperc
c762   format('Initial ',a8, ' configuration ',i3,' minimised D=', 
c     *      e13.6,' E%=',f6.2)      
      call encl(q)
       
       if(dopt.lt.doptm) then
       ntry=ntries
       doptm=dopt
	   dactm=dact
       initx=inital
      dy=ry(q+1)-ly(q+1)
	dx=rx(q+1)-lx(q+1)
	xsmall=lx(q+1)
	ysmall=ly(q+1)
      skinxx=0
c     write(*,*)'dx,dy=',dx,dy
	do  122 i=1,q
c
c establish  coordinates of rectangle corners in unit square.
c          
 
          x(1,i)=(lx(i)-xsmall)/dx
          y(1,i)=(ly(i)-ysmall)/dy
          x(2,i)=x(1,i)
          y(2,i)=(ry(i)-ysmall)/dy
          
          x(3,i)=(rx(i)-xsmall)/dx
          y(3,i)=y(2,i)
          x(4,i)=x(3,i)
          y(4,i)=y(1,i)
c  compute skinnyness   
      a=abs(lx(i)-rx(i))
      b=abs(ly(i)-ry(i))
      if (a .gt. b) then
       skin=(a/b)
      else
       skin=(b/a)
      end if

      if(skin.gt.skinxx) skinxx=skin
c      write(*,*)'skin,skinxx',skin,skinxx
 
122   continue
c  q+1 th is boundary 
c  rectangles numbered clockwise from lower left corner
c   i.e. 1=(0,0), 2(0,1), 3=(1,1), 4=(1,0)
          x(1,q+1)=0
          y(1,q+1)=0
          x(2,q+1)=0
          y(2,q+1)=1
          x(3,q+1)=1
          y(3,q+1)=1
          x(4,q+1)=1
          y(4,q+1)=0
      
      do 123 j=0,63
      xnopt(j)=xn(j)
      xmopt(j)=xm(j)
      ymopt(j)=ym(j)
123   continue
      
       end if  
c     endif of if(dopt.lt.doptm) then
c
c  try up to ntriex different starting positions
c  output the best in terms of optimality criterion 
      
        ntries=ntries + 1
      if( ntries.le.ntriex )then
      go to 100
      end if
c----------------------------------------------------------------------
 
101   continue
       call intpr(initx,-1,iarray,0)
c     write(*,658)ntry
c658   format('Best fit is from initial configuration ',i2)  
c       write(*,657) doptm
c657   format('Optimised criterion=', e13.6)
      
c   ** must delete intpr code for srd package **
c      write(string,564) initx
c      string =" Fitted from initial "//initx//" configuration"
c      call intpr(string,-1,iarray,0)
c      write(*,658) skinxx
c658   format('Maximum rectangle thinness=',f8.2)     
      
c  output  error
      err=0
      sumf=0
      suma=0
      ng=0
      xng=0
c&&      write(*,667)
      
c&&667   format('Cell Freq    Area    Freq-Area     Binary combination')  


      do 125 k=0,2**q-1 

         
       er=no(k) - xnopt(k)
c   ignore  whitespace i f no whitespace
       err=err+ abs(er)
c       
c   cell statistics to be returned to R  in c()
cccc       call realpr("no(k)",5,no(k),1)
        c(3,k+1)=no(k)
        c(4,k+1)=xnopt(k) 
        c(5,k+1)=100*zdiv(er,t)
        
c   cell coordinates to write labels         
        if(xnopt(k).gt.0) then
        c(1,k+1)=(xmopt(k) - xsmall)/dx
c  y position: slightly jitter? questionable benefit
c        c(2,k+1)=(ymopt(k) - ysmall)/dy  + 0.01*randu(iseed)
        c(2,k+1)=(ymopt(k) - ysmall)/dy 
        else 
c   no area

        c(1,k+1)=999
        c(2,k+1)=999
        end if
 
       
      call binno(k,q,bin)
      
c&&      if(.not.(int(no(k)).eq.0 .and. xnopt(k).le.t*0.0001) )then
c%%      write(*,666)k,
c&&     &      int(no(k)), xnopt(k), no(k)-xnopt(k), (bin(j),j=1,q)
c&&      end if
c&&666   format(i2, i8,2x,2f8.1,9x, 6i3)
      sumf=sumf+no(k)
      suma=suma+xnopt(k)
      
                do 126 i=1,q
           if(bin(i).eq.1)then
           ng(i)=ng(i)+no(k)   
           xng(i)=xng(i)+xn(k)
           end if
 
126     continue
      
      
125   continue
c     write(*,'(50(1h-))')
c     write(*,666) int(sumf), suma,sumf-suma
      e=100*(err/t)
c&&      write(*,*)'Error %=',e
c      if(e.gt.5) then
c     write(*,*)' Error>5%. Try re-running, may improve configuration'
c       end if
      
       
c      write(*,*)' Rectangle    Freq   Area'
       icheck=0
      do 1271 i=1,q
c     write(*,*) ng(i),xng(i)
      if(abs(ng(i)-xng(i)).gt.0.1)icheck=1
1271  continue
c      if(icheck.eq.1)write(*,*)'Note: Rectangle area-frequency unequal'
c     write(*,*)'Cell statistics:'
c      do k=0,2**q-1 
c     write(*,*) (c(j,k+1),j=1,5)
c          end do
c  return thin parameter (which will have been set in thin<0 on input) 
      
     
      thin=skinpa
      return
      end 
c  end of srd()      
c------------------------------------------------------------------      
      subroutine grandm(q,iseed)
c
c  establish randomised starting configuration
c   

      common/grpf/tot(6)
      Real tot
      integer q
      integer iseed
     

      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
    	common/space/whtesp,greysp,blcksp
c      common/tots/t,xt,tact
      
c  clock generated seed(-1) for random()      
c	call seed(-1) 
c	call seed(0) 
	

      s=sqrt(greysp+blcksp)
c      
c  position square rectangles inside unit square at random.
c  the square must be within a distance x of the edge of unit
c  square. But if subgroups are small, set a boundary within the
c  unit square. 

c  once boundary edge is fixed then considering X direction
c  there must be

c      |...x/2..|------------ s-x ------------|...x/2...|
c  for entire length to be s.
c  so left point of square must be uniform on s-x

      do 100 i=1,q
      x=sqrt(tot(i))

      side=(s-x)

c      call random(u)
      u=randu( iseed )     


	u=u-0.5 
	halfxx=u*side
      lx(i)=-halfxx
      
      rx(i)=lx(i) +x
c      call random(u)
       u=randu( iseed )     
	u=u-0.5 
	halfxx=u*side
      ly(i)=halfxx
	ry(i)=ly(i)- x
c     write(*,*)'start lx(i),ly(i),rx(i),ry(i)', lx(i),ly(i),rx(i),ry(i)
100    continue

      return
      end
c----------------------------------------------------------
       subroutine marg(q)
       common/grpf/t(6)
       Real t
       integer bin(6)
       integer q
       common/frqn/n(0:63)
       Real n
       do 3 j=1,q
3      t(j)=0
       do 1 i=1,2**q-1
       call binno(i,q,bin)
       do 2 j=1,q
2      if(bin(j).eq.1)t(j)=t(j)+n(i)
1      continue
c     write(*,*)'t=', t
       return
       end
c------------------------------------------------------
      double precision function func(p)
      double precision p(18)
      common/ngrp/q
      integer q
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      common/grpf/tot(6)

	double precision Dopt
      Real tot
      np=0
      do 10 l=1,q
      tq=tot(l)
      np=np+1
      lx(l)=p(np)
      np=np+1
      rx(l)=p(np)
      np=np+1
      ry(l)=p(np)
      ly(l)=ry(l)+zdiv(tq,rx(l)-lx(l))
10    continue
      call encl(q)
      call sectio(q,Dopt) 
      func=Dopt
      return
      end
cc---------------------------------

       subroutine encl(q)
c  enclose in unit square

      integer q  
      common/frqn/n(0:63)
      common/area/xn(0:63)
       Real xn,n
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
	common/inner/lxmin,lymax,rymin,rxmax
       Real lxmin,lymax
      Real lx,ly,rx,ry
	  logical centre
		centre= .TRUE. 

	vbig=1.0e10

      rxmax=-vbig
      rymin=vbig
      lxmin=vbig
      lymax=-vbig
      do 776 i=1,q
      if(rx(i).gt.rxmax)rxmax=rx(i)
      if(ry(i).lt.rymin)rymin=ry(i)
      if(lx(i).lt.lxmin)lxmin=lx(i)
      if(ly(i).gt.lymax)lymax=ly(i)
c  is possible the left and right switched?
c  
      if(lx(i).gt.rxmax)rxmax=lx(i)
      if(ly(i).lt.rymin)rymin=ly(i)
      if(rx(i).lt.lxmin)lxmin=rx(i)
      if(ry(i).gt.lymax)lymax=ry(i)
c
776   continue

c
c   align border, as q+1 th rectangle 
c 
c  NOTE: major (4 apr 2020) change here to force whtespace congruence
c  instead of total area congruence. tot is sum of 
c  "interior" areas, augmented by actual whitespace frequency
      tot=0
      do 144 j=1,2**q-1
       tot=tot+xn(j)
144    continue
c 
       tot = tot + n(0) 
	iq=q+1

c  new method to centre

	w=rxmax-lxmin
	d=lymax-rymin
	if(centre) then
c  following shouldnt happen      
	if(w.le.0.or.d.le.0) then 
	hw=sqrt(tot)
	hd=hw
	else
c
	hw=sqrt(tot*w/d)

	hd=tot/hw
	end if
	else
	hd=h
	hw=w
	end if
	

		hwside =hw*0.5 
		hdside= hd*0.5
	centx=(rxmax+lxmin)*0.5
	centy=(rymin+lymax)*0.5
c  force enclosure so that it  contains whole of rectangle configuration
	lx(iq)=min(centx-hwside,lxmin)
	rx(iq)=max(centx+hwside,rxmax)
      ry(iq)=min(centy-hdside,rymin)

c  place top and bottom border either against top
c  leaving space below

c  or place in middle with space above and below
c  by specifying yside*0.5
 
	 yside=zdiv(tot,(rx(iq)-lx(iq)))
c      ry(iq)=min(centy-yside*0.5,rymin)
      
	ly(iq)=max(yside+ry(iq),lymax)
c	ly(iq)=max(centrey+yside*0.5,lymax)  

       call reshap(q)


        return
        end
c  end encl
c-----------------------------------------------------

      subroutine reshap(q)
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      integer q

      i=q+1
      a=sqrt(abs(zdiv( (rx(i)-lx(i)),(ly(i)-ry(i)) )) )
      
      if(a.eq.0) a=1  
      do 777 i=1,q+1
      rx(i)=rx(i)/a
      lx(i)=lx(i)/a
      ry(i)=ry(i)*a
      ly(i)=ly(i)*a
777   continue
      return
      end 

c-----------------------------------------------------
       subroutine binno(n,q,bin)

       implicit integer(a-z)
       dimension bin(*)
c
c  to obtain binary representation of number n
c  in q 0/1 digits.
       pow=2**(q-1)
       rem=n
         do 1 j=q,1,-1
         i=rem/pow
         if(i.eq.1) rem=rem-pow
         bin(j)=i
         pow=pow/2
1        continue
       return
      end
c----------------------------------------------
c  inverse of binno - return n  for given bin() 
       subroutine ibinno(n,q,bin)
       implicit integer(a-z)
       dimension bin(*)

	n=0
	do 144 j=1,q
	if(bin(j).eq.1)n=n+2**(j-1)
144   continue
	return
      end
c---------------------------------------------------------
      subroutine wobble(q, funnew)
      integer q
      double precision p(18), xi(18,18)

      Real psave(18)
      common/grpf/tot(6)

      Real tot
      common/wobbl/ftol,itmax
	double precision ftol
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry

	double precision funnew, funold, func
      common/tots/t,xt,tact
       ftol=1.0e-10

 	itmax=20  
c  positive definite matrix for powell
c  18  = 6*3 max number of parameters q<=6  
c  I think  xi is same Order as par
      npx=18
      st=sqrt(t)
      do 1 i=1,18
      do 1 j=1,18
      xi(i,j)= 0.01*st 

      
      if(i.eq.j) xi(i,j)=0.1*st  
1     continue
c
c 
c     
      np=0
      do 10 l=1,q  

      np=np+1
      p(np)=lx(l)
      np=np+1
      p(np)=rx(l)
      np=np+1
      p(np)=ry(l)
10    continue
      do 11 i=1,np
      psave(i)=p(i)
11    continue

      funold=func(p)
	  
cccc	    call dblepr("funold",6,funold,1)
      call powell(p,xi,np,npx,ftol,iter,funnew,itmax)
c     write(*,754)funold,funnew, iter
c754   format(' Starting D=', e13.6,' Final D= ',e13.6,'  after ', i4,
c     * ' iterations')
cccc	    call dblepr("funnew",6,funnew,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      np=0
      do 107 l=1,q
      tq=tot(l)
      np=np+1
      lx(l)=p(np)
      np=np+1
      rx(l)=p(np)
      np=np+1
      ry(l)=p(np)
      ly(l)=ry(l)+zdiv(tq,rx(l)-lx(l))
107    continue


      return
      end  
c-----------------------------------------------
      
      subroutine sectio(q,Dopt)

      common/option/skinpa,absent, skinx

      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry


      common/object/lad,lsd,logl,chi2
      common/error/e
      logical lad,lsd,logl,chi2
      common/area/xn(0:63)
      common/frqn/n(0:63)
      common/posn/xm(0:63), ym(0:63)
      common/actual/dact
	  double precision dact
      Real n,xn
c   DIMENSION 14  =2*(6+1) q<=6
      Real x(14), y(14)
c    work array w() used in sort2 sorting routine
      integer w(14)
      integer q
      double precision xmid,ymid,area

      common/tots/t,xt,tact
	double precision Dopt
      Real areax(0:63)
      real widex(0:63)

      data small/1.0e-5/
c  evaluates rectangle area by partitioning into sections
c  with label=F return Dopt
c  with label=T  does not change xn array, but positions text
c  on diagram via outputlabels
	
      areax=0
      widex=0
      nds=2**q



      do 99 i=0,nds-1
      xn(i)=0.0
99    continue
c   form section arrays unsorted
      do 1 k=1,q+1
      k2=k*2
      x(k2-1)=lx(k)
      x(k2)=rx(k)
      y(k2-1)=ry(k)
      y(k2)=ly(k)
1     continue

      nn=2*(q+1)
      call sort2(nn,x,w)
	ncol=nn-1
      call sort2(nn,y,w)
	nrow=nn-1
c---------------------------------------------------------------
c eliminate cells outside border
	xlow=lx(q+1)
	xhigh=rx(q+1)
c	write(7,*)'xlow,xhigh',xlow,xhigh
	ncol=0
	do 144 l=1,nn
	if(x(l).ge.xlow) then 
	if(x(l).le.xhigh) then
	ncol=ncol+1
	x(ncol)=x(l)
	end if
	end if
144   continue
	ncol=ncol-1
c	write(7,*)'nn,ncol',nn, ncol

	ylow=ry(q+1)
	yhigh=ly(q+1)
c	write(7,*)'ylow,yhigh',ylow,yhigh
	nrow=0
	do 145 l=1,nn
	if(y(l).ge.ylow) then 
	if(y(l).le.yhigh) then
	nrow=nrow+1
	y(nrow)=y(l)
	end if
	end if
145   continue
	nrow=nrow-1
c	write(7,*)'nn,row',nn, row

c----------------------------------------------------------------
	nsubce=0

         do 2 k=1,ncol
c	write(7,*)'k,x(k),y(k)',k,x(k),y(k)
         xmid=(x(k+1)+x(k))*0.5
	   xd=x(k+1)-x(k)
         do 2 l=1,nrow
         ymid=(y(l+1)+y(l))*0.5
	   yd=y(l+1)-y(l)

         area=(x(k+1)-x(k))*(y(l+1)-y(l))
c         wide=abs(x(k+1)-x(k))

c	if( x(k+1)-x(k).le.0.05.or.y(l+1)-y(l).le.0.05) go to 2

c	if(label) write(7,*) 'Note empty grid cells'
c   ESTABLISHES the possibility of zero grid cell area
c	end if
c
c  determine which cell section s is in
         i=0
         do 3 j=1,q
c   midpoint
      if(xmid.lt.max(lx(j),rx(j)).and.xmid.gt.min(lx(j),rx(j)) )then
      if(ymid.lt.max(ly(j),ry(j)).and.ymid.gt.min(ly(j),ry(j))) then
      i=i+2**(j-1)
      
      end if
      end if
3        continue
c  keep coordinates where cell labels will be placed
c   by  largest sub-cell area
c or maybe by widest sub-cell? Later: area may be best         
         if(area.gt.areax(i)) then
             areax(i)=area
c         if(wide.gt.widex(i)) then
c          widex(i)=wide
         xm(i)=xmid
         ym(i)=ymid
         end if 
c   areas only inside a rectangle or bounding q+1 th 
c   rectangle allowed.
c--------------------------------------------------------------------   
c  contiguity of cells of type i?
c   required for contiguity routine. Only consider non-zero ares



c	if(label)write(7,*)'i ncells(i), k, l',i,ncells(i), k, l
c  try stripping off boundaries

c----------------------------------------------------------------

      if(i.gt.0.or.(i.eq.0.and.
     *    xmid.gt.lx(q+1).and.xmid.le.rx(q+1).and.
     *    ymid.gt.ry(q+1).and.ymid.le.ly(q+1)) )
     *    xn(i)=xn(i)+area

2       continue  
c  label 2  for k=1,ncol & l=1,row
c
c
c calculate objective function 
c	


	erx=0.
      Dopt=0.0
      xlogl=0
      t=0
       e=0

      do 244 i=0,nds-1
      t=t+n(i)

244   continue
c  arbitrary assign "absentcells" weight log(1/t)
      zeroc=-log(t)

      tx=0
      mssd=0
	small=t*0.0001
     
      do 450 i=0,nds-1

c  try a patch for possibility that enclosure is too small
c

	if(i.eq.0.and.xn(i).lt.n(i)) then
	diff=n(0)-xn(0)
	B=2*(abs(rx(q+1)-lx(q+1))+abs(ly(q+1)-ry(q+1)))
c additional boundary is soln to quadratic
c   4*delta**2 +B*delta=diff
	
c	delta=(-B+sqrt(B*B+4*4*diff))/(2*4)
	delta=(-B+sqrt(B*B+16*diff))/8
	rx(q+1)=rx(q+1)+delta
	lx(q+1)=lx(q+1)-delta


	ry(q+1)=ry(q+1)-delta
	ly(q+1)=ly(q+1)+delta
c	write(7,*)'patching bounadry by delta=',delta, n(0),xn(0)
	xn(i)=n(i)
	end if 

c  absent/missing cell
      if(xn(i).eq.0 .and. n(i).gt.0)mssd=mssd+1 
c  empty cell, maybe less important to avoid than missing??     
      if(xn(i).gt.small .and. n(i).eq.0)mssd=mssd+1 
      er= abs(xn(i)-n(i))
      e=e +er
	if(er.gt.erx) erx=er
      tx=tx+xn(i)


      if(lad) Dopt=Dopt+ er
      if(lsd) Dopt=Dopt+ er*er
      if(chi2) Dopt=Dopt+ er*er/max(0.5,xn(i))

	
	if(logl) then
         if(xn(i).gt.0) then
         xlogar=log(xn(i))
         else
         xlogar=zeroc
         end if 

      xlogl=xlogl+n(i)*xlogar
      end if 
      
450   continue

      if(logl) then
      if(tx.gt.0) then
c  need to scale by total area, to ensure areas are "proportions"
      xlogl=xlogl-t*log(tx)
      end if
      Dopt=-xlogl
      end if  
	  if(lsd) then
	  Dopt = sqrt(Dopt)
	  end if 
      e=100*e/t

c--------------------------------------
c  retain unpenalised statistic in common/actual/
      dact=Dopt
c penalty for absent and no data cells
      if(mssd.gt.0) then
c	  Dopt = Dopt *(nds/(nds-mssd))
		Dopt = Dopt *(mssd+1)
	  end if 

c  penalty for thinness
      pskin=1.0
      if(skinpa.gt.0.0) then
      skinx=0
      do 547 l=1,q
      a=abs(lx(l)-rx(l))
      b=abs(ly(l)-ry(l))
      if (a .gt. b) then
       skin=(a/b)
      else
       skin=(b/a)
      end if
      if(skin.gt.skinx) skinx=skin

547    continue
c	write(7,*) 'Dopt, skinx ',Dopt, t*skin

	  pskin=skinx**skinpa

c	  call realpr("pskin",5,pskin,1)
c	  call realpr("skinx",5,skinx,1)
      Dopt=Dopt + t*pskin 
c	write(7,*)'skinny penalising by', pskin
      end if !if(skinpa
      

      return
      end   
c-------------------------------------------------
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      Real arr(*)
      INTEGER brr(*),btemp
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      Real a,atemp
      if(n.le.0)return 
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        atemp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=atemp
        btemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=btemp
        if(arr(l+1).gt.arr(ir))then
          atemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l).gt.arr(ir))then
          atemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l+1).gt.arr(l))then
          atemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=btemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        atemp=arr(i)
        arr(i)=arr(j)
        arr(j)=atemp
        btemp=brr(i)
        brr(i)=brr(j)
        brr(j)=btemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
c        if(jstack.gt.NSTACK)write(*,*)'NSTACK too small in sort2$'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
c------------------------------------------------------
      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      double precision fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.e-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      double precision ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX)
      double precision xicom(NMAX),brent
      double precision f1dim
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
       xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END       
c-------------------------------------------------
      FUNCTION f1dim(x)
      INTEGER NMAX
      double precision f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      double precision pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END        
c-----------------------------------------------------------
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      double precision brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
c  passed as f1dim      
      EXTERNAL f  
      PARAMETER (ITMAX=50,CGOLD=.3819660,ZEPS=1.0e-5)
      INTEGER iter
      double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,r
      double precision tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      d=0.0
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
c**    pause      'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
c----------------------------------------------------------------
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-5)
      double precision dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
c    my addition
      if(ax.eq.bx) return  
      nits=0

1     nits=nits+1
      if(nits.ge.5) then
c  appears to be no minimum between entry points. assumne
c  constant. return with current inital values
      return
      end if
      if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        den=max(abs(q-r),TINY)
        if(q-r.lt.0) den=-den
c        u=bx-zdiv((bx-cx)*q-(bx-ax)*r,2.*den)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*den)
c***   u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        
c**        if((bx-u)*(u-cx).gt.0.)then
          if((bx.gt.u.and.u.gt.cx).or.(bx.lt.u.and.u.lt.cx))then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
c---------------------------------------------------
      SUBROUTINE  powell(p,xi,n,np,ftol,iter,fret,ITMAX)

      INTEGER iter,n,np,NMAX,ITMAX
      double precision fret,ftol,p(np),xi(np,np),func
      character*80 quest
      EXTERNAL func
      PARAMETER (NMAX=20)
CU    USES func,linmin
      INTEGER i,ibig,j
      double precision del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)

      fret=func(p)
      iter=0
      if(itmax.eq.0) then 
      return
      end if


       quest='Starting iterating'
c      write(*,*)quest


      do 11 j=1,n
        pt(j)=p(j)
11    continue
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      

      if(abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))then
c     write(quest,677) iter
c     write(7,*)'fp,fret,ftol', fp,fret,ftol

c  are changes in p() big
c??          do i=1,n
c??          if(xit(i).gt.1.0e-5) go to 688
c??          end do
c	write(quest,677) iter
c677  format('  Fitting converged in ', i2,' iterations')
c     write(*,*)quest

        return
        end if
c688   if(iter.eq.ITMAX) then
       if(iter.eq.ITMAX) then
c     write(quest,679) iter
c679   format(' No convergence after ',i3,' iterations')


      return
      end if
c     write(quest,678) iter, fret
c678   format(' Iteration=',i2,' Optimized D=',e13.6)

c     write(*,*)quest

      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END
c-----------------------------------------------------

      Real function zdiv(a,b)
      Real TINY, a, b
      tiny=1.0e-20
c  function to divide a/b and return 0 if b = 0
c	zdiv=0.0
      if(abs(b).gt.TINY) then
      zdiv=a/b
      else
      zdiv=0.0
      end if
      return
      end
c------------------------------------------------------
 	subroutine gnestd(q)
      common/grpf/tot(6)
c   nested rectangles, symmetrical
      Real tot
      integer q
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      xmax=0
      do 100 i=1,q
      x=sqrt(tot(i))
      if(x.gt.xmax) xmax=x
100   continue

      xmax2=xmax/2
      do 101 i=1,q
      x=sqrt(tot(i))
      halfx=x/2
c   centre them
c      lx(i)=-halfx
c      ly(i)=halfx
c      rx(i)=halfx
c      ry(i)=-halfx
c    left edge them
c  this done to ensure that nested rectangles
c   are not stuck when pertubated. This is same as 
c  as unused g5x routine. 
      lx(i)=-xmax2
      ly(i)=-xmax2 +x
      rx(i)=-xmax2 +x
      ry(i)=-xmax2
101   continue

      return
      end
c------------------------------------------------
	subroutine gdisjo(q)
      common/grpf/tot(6)
c  disjoint
      Real tot
      integer q
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
	common/space/whtesp,greysp,blcksp

	y=sqrt(tot(1)) 
      halfy=y/2
      lx(1)=-halfy
      ly(1)=halfy
      rx(1)=halfy
      ry(1)=-halfy

c   small gap between 
	delta=0
	xstart=halfy+delta
 
	x=tot(2)/y

      lx(2)=xstart
      ly(2)=halfy
      rx(2)=x+xstart
      ry(2)=-halfy

	
	alefto=whtesp
	do 100 i=3,q
	alefto=tot(i)+alefto
100   continue

	d=alefto/(rx(2)-lx(1))
	xstart=lx(1)
      do 101 i=3,q
	x=tot(i)/d

      lx(i)=xstart
      ly(i)=halfy +d
      rx(i)=x+xstart
      ry(i)=halfy
	xstart=x+xstart
101    continue

      ll=1
      lu=q
      return
      end

c-----------------------------------------------------------
       subroutine gstart(q,ntries,inital,iseed)
      integer q, ntries, iseed
      character*9 inital 
       if(ntries.eq.1) then
          inital=' nested'
          call gnestd(q)
          return
          end if
 
       if(ntries.eq.2) then
          inital=' disjoint'
          call gdisjo(q)
          return
       end if
       
      if(ntries.eq.3)then
          if(q.eq.3.or.q.eq.4) then
          if(q.eq.4) call g4i
          if(q.eq.3) call g3i
          inital='indep~ent'
          return
          else 
          call grandm(q,iseed)
          inital=' random'
           return
          end if
      end if
      
            if(ntries.eq.4)then
           if(q.eq.3.or.q.eq.4) then
          if(q.eq.4) call g4m2
          if(q.eq.3) call g3m2
          inital='estimated'
           return
           else 
          call grandm(q,iseed)
          inital=' random'
          return
           end if
            end if
            
                 if(ntries.ge.5) then
         call grandm(q,iseed)
          inital=' random'

                     return
                 end if
      return
      end
c------------------------------------------------------------
      subroutine g4m2
c  Method 2,  ( variations 4, 5 )
c  subroutine to set up rectangles for scaled subgroup diagrams
c  for 4 subgroups.
c  Method 2 uses the template in which all 2**4-1=15 cells are
c  represented and attempts to re-scale ensuring areas:
c  1234, 234, 123, 124, 134, 23 14 are correctly scaled and
c  hoping remainder are correct when rectangles comleted.
c
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      real lx,ly,rx,ry
      common/frqn/a0,a1,a2,a12,a3,a13,a23,a123,a4,a14,a24,a124
     *     ,a34,a134,a234,a1234,padd(48)
      common/grpf/tot(6)
      real tot
c
      t=(a1+a2+a3+a4+a12+a13+a14+a23+a24+a34
     *     +a123+a124+a234+a134+a1234)
      vmax=max(a123,a134)
      hmax=max(a234,a124)

c   intersection >1% of total
      if(a1234/t.gt.0.01) then 
      x=sqrt(a1234)
	
      xd=x
      else
	
c	ifail=1
	xd=1  
	x=1
c	return

        if(vmax.eq.0.and.hmax.eq.0) then
        x=0
        xd=0
        else
				if(vmax.gt.hmax) then
				 xd=sqrt(vmax)
					 if(a1234.gt.hmax) then
					 x=a1234/xd
					 else

					 x=sqrt(hmax)
					 end if
				  else
				  x=sqrt(hmax)
					 if(a1234.gt.vmax) then
					 xd=a1234/x
					 else
					 xd=sqrt(vmax)
					 end if
				   end if
		end if
	end if
         w=sqrt(a34)
         v=w
         u=sqrt(a12)
         z=u
 
         
      p=zdiv(a24,x)
      q=zdiv(a13,xd)

c total areas... (additional unrepresented areas a12 and a34 included)
      t1=tot(1)
      t2=tot(2)                                             
      t3=tot(3)
      t4=tot(4)
      r=max(0.0,zdiv(a14+a124,x+v)-z)
      s=max(0.0,zdiv(a23+a123,w+xd)-u)
      r1=zdiv(xd+z+r,u+x+v+q)
      w1=max(xd+z+r,sqrt(t1*r1))
c**      w1=xd+z
      d1=max(u+x+v,zdiv(t1,w1))
      w3=xd+w
      d3=s+u+x+v+q

      r2=zdiv(x+u+s,p+w+xd+z)
      d2=max(x+u+s,sqrt(t2*r2))
      w2=max(w+xd+z,zdiv(t2,d2))
      d4=x+v
      w4=p+w+xd+z+r
c  cordinated of left (l) and right (r) hand corners..
c wrt origin in centre of central intersection
      halfx=x/2
      halfxd=xd/2 
      lx(1)=-halfxd
      ly(1)=(d1-u-halfx)
      rx(1)=w1-halfxd
      ry(1)=-(halfx+u)
c
      lx(2)=-(w2-z-halfxd)
      ly(2)=halfx
      rx(2)=halfxd+z
      ry(2)=-d2+halfx
c
      lx(3)=-(w+halfxd)
      ly(3)=(q+v+halfx)
      rx(3)=halfxd
      ry(3)=-(halfx+u+s)
c
      lx(4)=-(halfxd+w+p)
      ly(4)=halfx+v
      rx(4)=halfxd +z+r
      ry(4)=-halfx
      return
      end
c-----------------------------------------------------
       subroutine g4i
c
c  create a q=4 diagram  assuming independence
c
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      common/grpf/tot(6)
      Real tot
      common/frqn/a0,a1,a2,a12,a3,a13,a23,a123,a4,a14,a24,a124
     *     ,a34,a134,a234,a1234,padd(48)
c
      t=(a0+a1+a2+a3+a4+a12+a13+a14+a23+a24+a34
     *     +a123+a124+a234+a134+a1234)
c
      t1=tot(1)
      t2=tot(2)                                             
      t3=tot(3)
      t4=tot(4)
c
       p1=t1/t
       p2=t2/t
       p3=t3/t
       p4=t4/t

	s=sqrt(t)
c
	lx(1)=0
	rx(1)=p1*s
	ly(1)=s
	ry(1)=0

	lx(2)=p1*(1-p2)*s
	rx(2)=lx(2)+p2*s
	ly(2)=s
	ry(2)=0

	lx(3)=0
	rx(3)=s
	ly(3)=p3*s
	ry(3)=0
c
	lx(4)=0
	rx(4)=s
	ry(4)=p3*(1-p4)*s
	ly(4)=ry(4)+p4*s
      return
      end
c-----------------------------------------------------
       subroutine g3i
c
c  create a q=3 diagram  assuming independence
c
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      common/grpf/tot(6)
      Real tot
      common/frqn/a0,a1,a2,a12,a3,a13,a23,a123,padd(56)
c
      t=(a0+a1+a2+a3+a12+a13+a23+a123)
c
      t1=tot(1)
      t2=tot(2)                                             
      t3=tot(3)
       p1=t1/t
       p2=t2/t
       p3=t3/t
c scale up to absolute size
      s=sqrt(t)
      

	lx(1)=0
	rx(1)=p1*s
	ly(1)=s
	ry(1)=0

	lx(2)=p1*(1-p2)*s
	rx(2)=lx(2)+p2*s
	ly(2)=s
	ry(2)=0

	lx(3)=0
	rx(3)=s
	ly(3)=p3*s
	ry(3)=0
c
      return
      end
c---------------------------------------------------------------    
         subroutine g3m2
c  
c  subroutine to set up rectangles for scaled subgroup diagrams
c  for 3 subgroups.
c
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      common/frqn/a0,a1,a2,a12,a3,a13,a23,a123,padd(56)
      common/grpf/tot(6)
      Real tot

c  total areas
      t1=tot(1)
      t2=tot(2)
      t3=tot(3)
      if(a123.gt.0.0) then
         if(a13.gt.0.and.max(a12,a23).gt.0) then
         s=sqrt(a13/(max(a12,a23)))
         else
         s=1
         end if
c
c   s=proportional scale of central intersection region (s=1 square)
c   fix intersection region as box x by xd
c
      x=sqrt(a123/s)
      xd=x*s
c
c  fix sizes of "wings" off the central region
c
      w=a23/x
      z=a12/x
      y=a13/xd
      else
c
c  use fix up for a123=0
c   
        xd=0
        w2=sqrt(t2) 
           if(a12.ne.0.or.a23.ne.0) then
           w=w2*(a23/(a12+a23))
           z=w2-w
           x=(a23+a12)/w2
           else
c
c   all three groups disjoint
c          a12=a23=a13=0
c               
           x=0
           z=w2/2
           w=w2/2
           end if
         y=max(0.0,sqrt(t3)-x)  
c      makes 3rd group  square ensures y>0
      end if    
c   end of if(a123=0)
c
c  set up width and depth of 3 rectangles
c
      w1=xd+z 
      d1=max(x+y, zdiv(t1,w1))
c
      w2=w+xd+z
      d2=max(x, zdiv(t2,w2))

      d3=x+y
      w3=max(w+xd, zdiv(t3,d3))
c
c  cordinated of left (l) and right (r) hand corners..
c wrt origin in centre of central intersection
      halfx=x/2
      halfxd=xd/2
      lx(1)=-halfxd
      ly(1)=d1-halfx
      rx(1)=z+halfxd 
      ry(1)=-halfx
c
      lx(2)=-w-halfxd
      ly(2)=halfx
      rx(2)=halfxd+z 
      ry(2)=-d2+halfx
c
      lx(3)=-w3+halfxd
      ly(3)=halfx+y
      rx(3)=halfxd
      ry(3)=-halfx
c
      return
      end  
c--------------------------------------------------
       subroutine g2(a1,a2,a12)
c  
c  subroutine to set up rectangles for scaled subgroup diagrams
c  for 2 subgroups, using offset construction.
c
      common/cnrs/lx(16),ly(16),rx(16),ry(16)
      Real lx,ly,rx,ry
      s=1.0
      if(a12.gt.0.0) then
c
c   s=proportional scale of central intersection region (s=1 square)
c   fix intersection region as box x by xd
c
      x=sqrt(a12/s)
      xd=x*s
c
      else
      x=0.0
      xd=0.0
      end if
c
c  solve quadratic for offset of area 1
      o1=(-x-xd+sqrt((x+xd)**2+4*a1) )/2.0
c  solve quadratic for offset of area 2
      o2=(-x-xd+sqrt((x+xd)**2+4*a2) )/2.0
c
c  cordinated of left (l) and right (r) hand corners..
c wrt origin in centre of central intersection
      halfx=x/2
      halfxd=xd/2
      lx(1)=-halfxd
      ly(1)=o1+halfx
      rx(1)=halfxd +o1
      ry(1)=-halfx
c
      lx(2)=-o2-halfxd
      ly(2)=halfx
      rx(2)=halfxd
      ry(2)=-o2-halfx
c
      return
      end
c----------------------------------------------------------
      function randu( seed )
c
cc randu returns a unit pseudorandom 0-1 uniform
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      randu= seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      randu
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, Real randu, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4huge
      parameter ( i4huge = 2147483647 )
      integer k
      integer seed
      Real randu

c      if ( seed .eq. 0 ) then
c        write ( *, '(a)' ) ' '
c        write ( *, '(a)' ) 'randu- Fatal error!'
c        write ( *, '(a)' ) '  Input value of SEED = 0.'
c        stop
c      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4huge
      end if

      randu= real( dble ( seed ) * 4.656612875D-10 )

      return
      end
