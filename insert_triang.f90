module precision

    implicit none
    integer, parameter :: DOUBLE = kind(1.0D+00)

end module precision


module constants

    use precision
    
    implicit none
    real(double), parameter :: ZERO  =  0.0D+00, &
                               ONE   =  1.0D+00, &
                               TWO   =  2.0D+00, &
                               THREE =  3.0D+00, &
                               FOUR  =  4.0D+00, &
                               FIVE  =  5.0D+00, &
                               SIX   =  6.0D+00, &
                               SEVEN =  6.0D+00, &
                               EIGHT =  8.0D+00, &
                               NINE  =  9.0D+00, &
                               TEN   = 10.0D+00, &
                               PI    = 3.1415926535897932D+00

end module constants


program main

    use precision
    use constants
    implicit real*8 (a-h,o-z)
!   ..    
!   .. Arguments ..    
    integer, allocatable :: elset(:), eltype(:), bcond(:,:), &
                            connec(:,:)

    real(double), allocatable :: coord(:,:)
!   ..    
!   .. Local arrays ..
    integer, allocatable :: holes(:,:), nnhole(:), elchk(:)
    integer :: lnods(4), polygon(20)
    character(len = 500) :: str
!   ..    
!   .. Parameters ..
    integer, parameter :: inpf = 10, outf = 20
    integer, parameter :: nnodhole = 20
!    
!   ..    
!       
    open(unit = inpf, file = 'inp_cb.dat', status = 'unknown')
    open(unit = outf, file = 'out_cb.dat', status = 'unknown')

    read(inpf,*) nnode, nelem, ndime

    call alloc_memory( nelem, nnode, ndime, &
         elset, connec, bcond, eltype, coord )

    allocate( holes(nnodhole,nelem) )
    allocate( nnhole(10*nelem) )
    allocate( elchk(nelem) )
    
!=====================================================================
!   Read from input file
!=====================================================================
    do inode = 1, nnode
        read(inpf,*) ii, (coord(jj,inode), jj = 1, ndime), (bcond(kk,inode), kk = 1, 2)
    end do

    nnope = ndime + 1
    do ielem = 1, nelem
        read(inpf,*) ii, elset(ielem), eltype(ielem), (connec(jj,ielem), jj = 1, nnope)
    end do

!=====================================================================
!   New functionalities
!=====================================================================

    ihole = 0
    do ielem = 1, nelem
        
        if ( elset(ielem) >= 4 .and. &
             elset(ielem) <= 5 .and. &
             elchk(ielem) == 0 ) then
            
!---------------------------------------------------------------------
            ihole = ihole + 1
            elchk(ielem) = 1
            
            lnods(1:3) = connec(1:3,ielem)
            lnods(4)   = connec(1,ielem)
!
!           Compute maximum distance between nodes in a HAR-element
            max_dist = ZERO
            do inode = 1, nnope
                dist = comp_dist( coord, lnods(inode), lnods(inode+1) )
                if ( dist > amax_dist ) amax_dist = dist
            end do

            ii = 2
            do inode = 1, nnope
                dist = comp_dist( coord, lnods(inode), lnods(inode+1) )
                if ( dist < 0.1 * amax_dist ) then

                     holes(ii-1,ihole) = lnods(inode)
                     holes(ii,ihole) = lnods(inode+1)
                     ! write(*,*) holes(ii-1,ihole)
                     ! write(*,*) holes(ii,ihole)
                     ! write(*,*) 
                     exit
                     
                 end if

            end do
!---------------------------------------------------------------------

            node0 = lnods(inode)
            
100         continue
            knode = holes(ii,ihole)
            do jelem = 1, nelem

                if ( elset(jelem) >= 4 .and. &
                     elset(jelem) <= 5 .and. &
                     elchk(jelem) == 0 ) then

                    lnods(1:3) = connec(1:3,jelem)
                    lnods(4)   = connec(1,jelem)

                    amax_dist = ZERO
                    do jnode = 1, nnope
                        dist = comp_dist( coord, lnods(jnode), lnods(jnode+1) )
                        if ( dist > amax_dist ) amax_dist = dist
                    end do
                    
                    do jnode = 1, nnope

                        dist = comp_dist( coord, lnods(jnode), lnods(jnode+1) )
                        if ( dist < 0.2 * amax_dist .and. &
                             ( lnods(jnode) == knode .or. &
                               lnods(jnode) == node0  ) &
                            ) then
                            
                            elchk(jelem) = 1

                            ii = ii + 1
                                    
                            holes(ii,ihole) = lnods(inode)
                            if ( lnods(jnode) == knode ) then

                                if ( lnods(inode) == knode ) &
                                    holes(ii,ihole) = lnods(inode+1)
                                
                                if (holes(ii,ihole) == node0) then
                                    nnhole(ihole) = ii
                                    go to 200
                                end if
                                                                
                                go to 100
                                
                            end if
                                        
                        end if
                            
                    end do
                    
                end if

            end do
            
!---------------------------------------------------------------------
            
        end if
200     continue
        
    end do

300 continue
    new_nnode = nnode
    new_nelem = nelem

    do ii = 1, ihole

!        write(*,*) nnhole(ii)
        polygon(1:nnhole(ii)) = holes(1:nnhole(ii),ii)
!        write(*,*) ( polygon(jnode), jnode = 1, nnhole(ii) ) 
        
        new_nnode = new_nnode + 1
        call centroid( (nnhole(ii)-1), coord, polygon, xc, yc )
        
        coord(1,new_nnode) = xc
        coord(2,new_nnode) = yc
        
        do ielem = 1, nnhole(ii) - 1  
            new_nelem = new_nelem + 1

            inode = ielem
            jnode = ielem+1
            
            connec(1,new_nelem) = new_nnode 
            connec(2,new_nelem) = polygon(inode) 
            connec(3,new_nelem) = polygon(jnode) 

            ! write(*,*) connec(1:3,new_nelem) 
            ! write(*,*) 
            
            lnods(1:3) = connec(1:3,new_nelem)
            lnods(4) = connec(1,new_nelem)
            
            area = area_polygon( 3, coord, lnods )

            if ( area < ZERO ) then
                iswap = connec(2,new_nelem)
                connec(2,new_nelem) = connec(3,new_nelem)
                connec(3,new_nelem) = iswap
            end if
            
            elset(new_nelem) = 6
            eltype(new_nelem) = 1

        end do
    end do
    
    ! do ii = 1, ihole
!        write(*,*) ( polygon(jnode), jnode = 1, nnhole(ii) ) 
    ! end do
        
!=====================================================================
!   Write to output file
!=====================================================================
    do inode = 1, new_nnode
        write(outf,"(i10,2(x,e20.13),2i10)") inode, (coord(jj,inode), jj = 1, ndime), (bcond(kk,inode), kk = 1, 2)
    end do

    do ielem = 1, new_nelem
        write(outf,"(6i10)") ielem, elset(ielem), eltype(ielem), (connec(jj,ielem), jj = 1, ndime + 1)
    end do

    close(inpf)
    close(outf)

!=====================================================================
!   Write to Paraview file
!=====================================================================

    open(unit = outf, file = 'para.vtu', status = 'unknown')
    write(outf,500) new_nnode, new_nelem
500 format('<?xml version="1.0"?>'/&
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'/&
        '<UnstructuredGrid>'/&
        '<Piece NumberOfPoints="', i10,'" NumberOfCells="', i10,'">'/&
        '<Points>'/&
        '<DataArray type="Float64"  NumberOfComponents="3" format="ascii">')
    
    do ino = 1, new_nnode
        write(outf,'(3(" ",e20.13))') &
            (coord(jj,ino), jj = 1, 2), 0.0D+00
    end do
    write(outf,510)
510 format('</DataArray>',/,'</Points>',/,'<Cells>')

    
520 format('<DataArray type="Int64"       Name="connectivity" format="ascii">')
    write(outf,520)
    do iel = 1, new_nelem

        nnel = 3        
        do ino = 1, nnel
            lnods(ino) = connec(ino,iel) - 1
        end do
        
        write(str,'("(",i5,"i10)")') nnel
        write(outf,str) (lnods(jj), jj = 1, nnel)        

    end do
    write(outf,530)
530 format('</DataArray>')

540 format('<DataArray type="Int64"            Name="offsets" format="ascii">')
    write(outf,540)
    kk = 0
    do iel = 1, new_nelem

        nnel = 3
        kk = kk + nnel
        write(outf,*) kk

    end do
    write(outf,530)
    
550 format('<DataArray type="Int64"            Name="types" format="ascii">')
    write(outf,550)
    do iel = 1, new_nelem
        lty = 5
        write(outf,'(i5)') lty 
    end do
    write(outf,530)

    write(outf,560)
560 format('</Cells>',/&
        '<CellData>',/&
        '<DataArray type="Int64" Name="Material" format="ascii">')
    do iel = 1, new_nelem
        write(outf,'(i5)') elset(iel)
    end do
    write(outf,530)
    
    write(outf,570)
570 format('</CellData>',/&
        '</Piece>',/&
        '</UnstructuredGrid>',/&
        '</VTKFile>')
    write(*,'(A)') ""
    close(outf)

    deallocate(holes)
!    deallocate(nnhole)
    deallocate(elchk)
    !deallocate( elset )
    deallocate( connec )
    deallocate( bcond )
    !deallocate( eltype )
     deallocate( coord )
!    call free_memory( elset, connec, bcond, eltype, coord )

    
contains


    subroutine alloc_memory( nelem, nnode, ndime, &
               elset, connec, bcond, eltype, coord )

        implicit none
        integer :: nelem, nnode, ndime
        integer, allocatable :: elset(:), eltype(:), bcond(:,:), connec(:,:)
        real(DOUBLE), allocatable :: coord(:,:)
        
        allocate(coord(ndime,2*nnode))
        allocate(bcond(2,2*nnode))

        allocate(elset(10*nelem))
        allocate(connec(4,10*nelem))
        allocate(eltype(10*nelem))
        
    end subroutine alloc_memory


    subroutine free_memory( elset, connec, bcond, eltype, coord )

        implicit none
        integer, allocatable :: elset(:), eltype(:), bcond(:,:), connec(:,:)
        real(DOUBLE), allocatable :: coord(:,:)
        
        deallocate( elset )
        deallocate( connec )
        deallocate( bcond )
        deallocate( eltype )
        deallocate( coord )
        
    end subroutine free_memory


    function comp_dist(coord, node1, node2)
        real(kind(1.0D+00)), intent(in) :: coord(:,:)
        integer, intent(in) :: node1, node2
        real(kind(1.0D+00)) :: comp_dist

        comp_dist = sqrt( &
            ( coord(1,node2) - coord(1,node1) )**2 + &
            ( coord(2,node2) - coord(2,node1) )**2 &
            )
        
        return
    end function comp_dist

    
    function area_polygon( nnel, coord, lnod ) result ( area )

        implicit none
        integer, intent(in) :: nnel
        integer, intent(in) :: lnod(nnel+1)
        real(DOUBLE), intent(in), allocatable :: coord(:,:)

        integer :: inode
        integer :: ii, jj
        real(DOUBLE) :: x1, y1, x2, y2, area
        
        area = ZERO
        do inode = 1, nnel

            ii = lnod(inode)
            jj = lnod(inode+1)
            
            x1 = coord(1,ii)
            y1 = coord(2,ii)

            x2 = coord(1,jj)
            y2 = coord(2,jj)

            area = area + (x1*y2 - y1*x2)

        end do
        area = area/TWO
        
    end function area_polygon
    


    subroutine centroid( nnel, coord, lnod, xc, yc )

        implicit none
        integer, intent(in) :: nnel
        integer, intent(in) :: lnod(nnel+1)
        real(DOUBLE), intent(in), allocatable :: coord(:,:)
        real(DOUBLE), intent(out) :: xc, yc

        integer :: inode
        integer :: ii, jj
        real(DOUBLE) :: x1, y1, x2, y2, area
        
        xc = ZERO
        yc = ZERO
        do inode = 1, nnel

            ii = lnod(inode)
            jj = lnod(inode+1)
            
            x1 = coord(1,ii)
            y1 = coord(2,ii)

            x2 = coord(1,jj)
            y2 = coord(2,jj)

            xc = xc + x1
            yc = yc + y1
            
            ! xc = xc + (x1 + x2)*(x1*y2 - x2*y1)
            ! yc = yc + (y1 + y2)*(x1*y2 - x2*y1)

         ! write(*,*) x1, y1
        end do
         ! write(*,*)             
        ! xc = xc + x2
        ! yc = yc + x2

        xc = xc/nnel
        yc = yc/nnel
        
        ! area = area_polygon( nnel, coord, lnod )
        ! xc = xc/(SIX*area)
        ! yc = yc/(SIX*area)

        ! write(*,*)
        ! write(*,*) xc, yc
        ! stop
        
    end subroutine centroid

    
end program main


