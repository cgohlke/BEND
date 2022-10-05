!Calculation of Bending and Curvature
!for Arbitrary DNA Base Sequences
!
!Version 2.0, January 1995
!
!David S. Goodsell
!Assistant Member
!The Scripps Research Institute
!goodsell@scripps.edu
!
!please reference the following report in any publication resulting
!from use of this program.
!Goodsell, D.S. & Dickerson, R.E. (1994) "Bending and Curvature Calculations
!in B-DNA" Nucl. Acids Res. 22, 5497-5503
!
!Input file: 'sequence.dat'
!            Ten bases per line
!            Maximum of 500 bases - ten lines
!example:
!CGCTCGTGCG
!CGCTCGTGCG
!CGCAAAAAAA
!CGCAAAAAAA
! etc.
!Input file: 'angles.dat'
!            Twist, roll and tilt angles for each step
!            file contains three comment lines, followed by
!            three 4X4 matrices
!            last line: rbend, rcurve, bendscale, curvescale
!            rbend sets the window for calculating bend values
!              the bend at base=n is the difference in base normals at
!              base=n+rbend and base=n-rbend
!              rbend should be less than 5
!            rcurve sets the window for calculating curvature similarly
!              rcurve should be greater than 10
!            bendscale and curvescale are multiplicitive scales for
!              drawing the profiles - I use 2.1 for bendscale and
!              0.7 for curvescale, i.e. bends scaled by three relative
!              to the curves
!example:
!**non A-tract bending model**
!Data from Calladine, Drew and McCall (1988) JMB 201, 127
! A   T   G   C  <--second base
!35.,34.,34.,34.  A first base ----------------------------
!34.,35.,34.,34.  T   |                              TWIST
!34.,34.,34.,34.  G   V
!34.,34.,34.,34.  C ---------------------------------------
!0.0,3.3,3.3,3.3  A ---------------------------------------
!6.6,0.0,3.3,3.3  T  (<-- 6.6 is roll for a TpA step, for example)
!3.3,3.3,3.3,3.3  G                                  ROLL
!3.3,3.3,3.3,3.3  C ---------------------------------------
!0.0,0.0,0.0,0.0  A ---------------------------------------
!0.0,0.0,0.0,0.0  T                                  TILT
!0.0,0.0,0.0,0.0  G
!0.0,0.0,0.0,0.0  C ---------------------------------------
!1.0,15.0,2.1,0.7
!
!Output file: 'profile.ps'
!             postscript file with sequence,
!             Bend profile in solid line,
!             Curve profile in dotted line
!
!Output file: 'profile.log'
!             ascii file with sequence, bend angle (degrees),
!             and curvature angle (degrees)
!
!-----------------------------------------------------------
		character sequence(500)
		integer isequence(500)
		real*4 twist(4,4),roll(4,4),tilt(4,4)
		real*4 dx,dy,x(500),y(500),pi,sc
		real*4 xave(500),yave(500)
		real*4 bend(500),bendscale,curve(500),curvescale

		sc=3.
		pi=3.14159

		do i=1,500
		sequence(i)=' '
		enddo

		open(1,file='sequence.dat',form='FORMATTED',status='OLD')
		open(2,file='angles.dat',form='FORMATTED',status='OLD')
		open(3,file='profile.log',form='FORMATTED',
     &		carriagecontrol='LIST',status='NEW')
		open(4,file='profile.ps',form='FORMATTED',
     &		carriagecontrol='LIST',status='NEW')

!Read twist, roll and tilt angles for each step
		read(2,*)
		read(2,*)
		read(2,*)
		do i=1,4
		read(2,*) (twist(i,j),j=1,4)
		do j=1,4
		twist(i,j)=twist(i,j)/180.*3.14159
		enddo
		enddo
		do i=1,4
		read(2,*) (roll(i,j),j=1,4)
		enddo
		do i=1,4
		read(2,*) (tilt(i,j),j=1,4)
		enddo
		read(2,*) rbend,rcurve,bendscale,curvescale
		if (rbend.gt.5) then
		  write(6,*) 'rbend>5, setting to 5'
		  rbend=5.
		endif

!Read in nucleotide sequence
		n=1
 10		read(1,11,end=99) (sequence(i),i=n,n+9)
 		write(6,100) (sequence(i),i=n,n+9)
 		n=n+10
 		goto 10
 11		format(10a1)
 100	format(2x,10a1)

 99		n=n-1
 		if (n.gt.500) then
 		  write(6,*) ' too many bases in input file : ',n
 		  write(6,*) ' truncate to 500 bases '
 		  stop
 		endif
 		do i=1,n
 		isequence(i)=4
 		if (sequence(i).eq.'A') isequence(i)=1
 		if (sequence(i).eq.'T') isequence(i)=2
 		if (sequence(i).eq.'G') isequence(i)=3
 		if (sequence(i).eq.'C') isequence(i)=4
 		enddo

		write(6,*) 'number of bases : ',n

!Calculate trajectory of helix axis
 		do i=1,n-1
 		twistsum=twistsum+twist(isequence(i),isequence(i+1))
 		dx=roll(isequence(i),isequence(i+1))*sin(twistsum)+
     &     tilt(isequence(i),isequence(i+1))*sin(twistsum-pi/2.) 
 		dy=roll(isequence(i),isequence(i+1))*cos(twistsum)+
     &     tilt(isequence(i),isequence(i+1))*cos(twistsum-pi/2.) 
 		x(i+1)=x(i)+dx
 		y(i+1)=y(i)+dy
 		enddo

!Calculate average of trajectory over 10 base pairs
 		do i=6,n-5
 		 rxsum=0
 		 rysum=0
 		 do j=-4,4
 		  rxsum=rxsum+x(i+j)
 		  rysum=rysum+y(i+j)
 		 enddo
 		  rxsum=rxsum+x(i+5)/2.
 		  rysum=rysum+y(i+5)/2.
 		  rxsum=rxsum+x(i-5)/2.
 		  rysum=rysum+y(i-5)/2.

 		  xave(i)=rxsum/10.
 		  yave(i)=rysum/10.
 		enddo

!Calculate bend profile
 		do i=rbend+1,n-rbend
 		  bend(i)=sqrt((x(i+rbend)-x(i-rbend))**2+
     & 		            (y(i+rbend)-y(i-rbend))**2)*bendscale

 		enddo

!Calculate curve profile
 		do i=rcurve+6,n-rcurve-5
 		curve(i)=sqrt((xave(i+rcurve)-xave(i-rcurve))**2+
     & 		           (yave(i+rcurve)-yave(i-rcurve))**2)
		enddo

		do i=1,n,50
		write(6,200) (int(curve(j)/10),j=i,i+49)
		write(6,200) (int(curve(j)-10*
     &                int(curve(j)/10)),j=i,i+49)
		write(6,201) (sequence(j),j=i,i+49)
		write(6,*)
		enddo

		do i=1,n
		write(3,202) sequence(i),bend(i),curve(i)
		enddo

 		do i=rcurve+6,n-rcurve-5
 		curve(i)=curve(i)*curvescale
		enddo

 200	format(1x,50i1)
 201	format(1x,50a1)
 202	format(1x,a1,3x,f6.1,3x,f6.1)

!Create profile in postscript
		write(4,1) '%!'
		write(4,*) ' /Courier findfont 15 scalefont setfont'
		write(4,*) ' 0.0 setgray'
 1		format(a2)

!Draw curve profile
		icurve=int(rcurve)
		xstart=104.
		rstep=9.
		yline=790.
		do i=1,n,50
		  yline=yline-70.
		  write(4,*) ' 0.5 setlinewidth'
		  write(4,*) ' [] 0 setdash'
		  write(4,502) xstart-4,yline,'moveto'
		  write(4,502) xstart+446,yline,'lineto'
		  write(4,*) ' stroke'
		  write(4,*) ' 1.5 setlinewidth'
		  write(4,*) ' [5 4] -2 setdash'
		  write(4,502) xstart,yline+curve(i),'moveto'
		  do j=0,49
		    r=float(j)
		    if ((i+j.gt.icurve+5).and.(i+j.lt.n-icurve-5)) then
		    if (i+j.gt.icurve+6) then
		     write(4,502) xstart+r*rstep,
     &	         yline+curve(i+j),'lineto'
     		else
		     write(4,502) xstart+r*rstep,
     &	         yline+curve(i+j),'moveto'
		    endif
		    endif
		  enddo
		  write(4,*) ' stroke'
		enddo

!Write base sequence
		xstart=100.
		yline=780.
		do i=1,n,50
		yline=yline-70.
		write(4,*) ' [] 0 setdash'
		write(4,502) xstart,yline,'moveto'
		 write(4,500) '(',(sequence(j),j=i,i+49),') show'
		enddo

!Draw bend profile
		write(4,*) ' 1.0 setlinewidth'
		xstart=104.
		yline=790.
		do i=1,n,50
		  yline=yline-70.
		  write(4,502) xstart,yline+bend(i),'moveto'
		  do j=0,49
		    r=float(j)
		    if ((sequence(i+j).ne.' ').and.(i+j.le.n-rbend)) then
		    if (i+j.gt.rbend+1) then
		     write(4,502) xstart+r*rstep,
     &	         yline+bend(i+j),'lineto'
     		else
		     write(4,502) xstart+r*rstep,
     &	         yline+bend(i+j),'moveto'
		    endif
		    endif
		  enddo
		  write(4,*) ' stroke'
		enddo

		write(4,501) 'showpage'
 500	format(1x,a1,50a1,a6)
 501	format(1x,a8)

 502	format(1x,2f8.2,1x,a6)
 503	format(1x,2f8.2,1x,a19)

 		end
