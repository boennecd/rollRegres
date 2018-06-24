      subroutine dchdd(r,ldr,p,x,z,ldz,nz,y,rho,c,s,info)
      integer ldr,p,ldz,nz,info
      double precision r(ldr,p),x(p),z(ldz,nz),y(nz),s(p)
      double precision rho(nz),c(p)
c
c     By B. Christoffersen. I downloaded the files from www.netlib.org/linpack.
c     for some reason all the dimension where 1 instead of nz or p.
c
c     dchdd downdates an augmented cholesky decomposition or the
c     triangular factor of an augmented qr decomposition.
c     specifically, given an upper triangular matrix r of order p,  a
c     row vector x, a column vector z, and a scalar y, dchdd
c     determineds a orthogonal matrix u and a scalar zeta such that
c
c                        (r   z )     (rr  zz)
c                    u * (      )  =  (      ) ,
c                        (0 zeta)     ( x   y)
c
c     where rr is upper triangular.  if r and z have been obtained
c     from the factorization of a least squares problem, then
c     rr and zz are the factors corresponding to the problem
c     with the observation (x,y) removed.  in this case, if rho
c     is the norm of the residual vector, then the norm of
c     the residual vector of the downdated problem is
c     dsqrt(rho**2 - zeta**2). dchdd will simultaneously downdate
c     several triplets (z,y,rho) along with r.
c     for a less terse description of what dchdd does and how
c     it may be applied, see the linpack guide.
c
c     the matrix u is determined as the product u(1)*...*u(p)
c     where u(i) is a rotation in the (p+1,i)-plane of the
c     form
c
c                       ( c(i)     -s(i)     )
c                       (                    ) .
c                       ( s(i)       c(i)    )
c
c     the rotations are chosen so that c(i) is double precision.
c
c     the user is warned that a given downdating problem may
c     be impossible to accomplish or may produce
c     inaccurate results.  for example, this can happen
c     if x is near a vector whose removal will reduce the
c     rank of r.  beware.
c
c     on entry
c
c         r      double precision(ldr,p), where ldr .ge. p.
c                r contains the upper triangular matrix
c                that is to be downdated.  the part of  r
c                below the diagonal is not referenced.
c
c         ldr    integer.
c                ldr is the leading dimension fo the array r.
c
c         p      integer.
c                p is the order of the matrix r.
c
c         x      double precision(p).
c                x contains the row vector that is to
c                be removed from r.  x is not altered by dchdd.
c
c         z      double precision(ldz,nz), where ldz .ge. p.
c                z is an array of nz p-vectors which
c                are to be downdated along with r.
c
c         ldz    integer.
c                ldz is the leading dimension of the array z.
c
c         nz     integer.
c                nz is the number of vectors to be downdated
c                nz may be zero, in which case z, y, and rho
c                are not referenced.
c
c         y      double precision(nz).
c                y contains the scalars for the downdating
c                of the vectors z.  y is not altered by dchdd.
c
c         rho    double precision(nz).
c                rho contains the norms of the residual
c                vectors that are to be downdated.
c
c     on return
c
c         r
c         z      contain the downdated quantities.
c         rho
c
c         c      double precision(p).
c                c contains the cosines of the transforming
c                rotations.
c
c         s      double precision(p).
c                s contains the sines of the transforming
c                rotations.
c
c         info   integer.
c                info is set as follows.
c
c                   info = 0  if the entire downdating
c                             was successful.
c
c                   info =-1  if r could not be downdated.
c                             in this case, all quantities
c                             are left unaltered.
c
c                   info = 1  if some rho could not be
c                             downdated.  the offending rhos are
c                             set to -1.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dchdd uses the following functions and subprograms.
c
c     fortran dabs
c     blas ddot, dnrm2
c
      integer i,ii,j
      double precision a,alpha,azeta,norm,dnrm2
      double precision ddot,t,zeta,b,xx
c
c     solve the system trans(r)*a = x, placing the result
c     in the array s.
c
      info = 0
      s(1) = x(1)/r(1,1)
      if (p .lt. 2) go to 20
      do 10 j = 2, p
         s(j) = x(j) - ddot(j-1,r(1,j),1,s,1)
         s(j) = s(j)/r(j,j)
   10 continue
   20 continue
      norm = dnrm2(p,s,1)
      if (norm .lt. 1.0d0) go to 30
         info = -1
      go to 120
   30 continue
         alpha = dsqrt(1.0d0-norm**2)
c
c        determine the transformations.
c
         do 40 ii = 1, p
            i = p - ii + 1
            scale = alpha + dabs(s(i))
            a = alpha/scale
            b = s(i)/scale
            norm = dsqrt(a**2+b**2+0.0d0**2)
            c(i) = a/norm
            s(i) = b/norm
            alpha = scale*norm
   40    continue
c
c        apply the transformations to r.
c
         do 60 j = 1, p
            xx = 0.0d0
            do 50 ii = 1, j
               i = j - ii + 1
               t = c(i)*xx + s(i)*r(i,j)
               r(i,j) = c(i)*r(i,j) - s(i)*xx
               xx = t
   50       continue
   60    continue
c
c        if required, downdate z and rho.
c
         if (nz .lt. 1) go to 110
         do 100 j = 1, nz
            zeta = y(j)
            do 70 i = 1, p
               z(i,j) = (z(i,j) - s(i)*zeta)/c(i)
               zeta = c(i)*zeta - s(i)*z(i,j)
   70       continue
            azeta = dabs(zeta)
            if (azeta .le. rho(j)) go to 80
               info = 1
               rho(j) = -1.0d0
            go to 90
   80       continue
               rho(j) = rho(j)*dsqrt(1.0d0-(azeta/rho(j))**2)
   90       continue
  100    continue
  110    continue
  120 continue
      return
      end
