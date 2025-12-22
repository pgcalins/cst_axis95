!
!     cst_axis95  Elementos finitos elasticos axissimetricos cst =
!     ==========  triangulo de deformacao constante
!
!     Adaptado de:
!     BREBBIA,C.A. & FERRANTE,A.J. (1986) Computational methods for the
!     solution of engineering problems. London: Pentech Press. 370p.
!
!     BERGERON, William Joseph (1968) Finite Element Analysis of Salt
!     Pillar Models.  LSU Historical Dissertations and Theses. 1466.
!     https://digitalcommons.lsu.edu/gradschool_disstheses/1466
!
!     P¢s-processamento no GiD https://www.gidhome.com/
!
!     Implementado por:
!     Paulo Gustavo Cavalcante Lins <pgcalins@gmail.com>
!
module variaveis
implicit none

  integer, parameter :: ndf=2 !! numero de graus de liberdade por n¢
  integer, parameter :: nne=3 !! numero de n¢s por elemento
  integer, parameter :: ndfel=ndf*nne !! graus de liberdade por elemento

  integer, parameter :: in=15 !! numero do arquivo de entrada
  integer, parameter :: io=16 !! numero do arquivo de saida
  integer, parameter :: i33=33 !! numero do arquivo do GiD
  integer, parameter :: i34=34 !! numero do arquivo do GiD

  integer :: nnode !! numero de n¢s
  real(4), allocatable :: X(:),Y(:) !! coordenadas dos n¢s
  integer :: nelem !! numero de elementos
  integer, allocatable :: kon(:) !! conectividade dos elementos
  integer, allocatable :: imat(:) !! numero do material do elemento
  integer :: nmat !! numero de materiais
  real(4), allocatable :: Ei(:) !! modulo de elasticidade do elemento
  real(4), allocatable :: Poisson(:) !! Poisson do elemento
  integer :: nln !! numero de nos carregados
  integer :: nbn !! numero de nos com condicao de contorno
  integer, allocatable :: istatus(:) !! indicador de status
  real(4), allocatable :: Prescrito(:) !! deslocamentos prescritos

  integer :: neq !! numero total de incognitas
  real(4), allocatable :: Carga(:) !! Vetor de cargas nodais
  real(4), allocatable :: Desloc(:) !! Vetor de deslocamentos nodais
  real(4), allocatable :: RigGlobal(:,:) !! Matriz de rigidez global
  real(4) :: RigElem(ndfel,ndfel) !! Matriz de rigidez do elemento
  real(4) :: DB(4,6) !! Matriz da lei de Hooke multiplicada por B
  real(4) :: ORX,ORY !! Centroide do elemento corrente

  real(4), allocatable :: stress(:,:) !! Tensoes no centroide dos elementos

end module variaveis



module Entrada_Saida

contains

subroutine Abre_Arquivos()
use variaveis

   open(UNIT=in,file='NomeArq.dat',status='old')
   open(UNIT=io,file='NomeArq.out',status='UNKNOWN')
   open(UNIT=i33,file='NomeArq.post.msh',status='UNKNOWN')
   open(UNIT=i34,file='NomeArq.post.res',status='UNKNOWN')

   return
end subroutine Abre_Arquivos



subroutine Entrada_Dados()
use variaveis
implicit none
integer :: i,j,ic(nne),n1,k,k1,L1,L2,n2
real(4) :: w(ndf)

   write(io,'(A,/)')' DADOS DE ENTRADA'

   read(in,*)nnode
   write(io,'(A,i5)')' numeros de nos          :',nnode
   write(io,'(/,A)')' coordenadas nodais'
   write(io,'(7x,A,6x,A,9x,A)')' no ','x','y'
   allocate(X(nnode)); allocate(Y(nnode))
   do j=1,nnode
      read(in,*)i,X(i),Y(i)
      write(io,'(i10,2f10.2)')i,X(i),Y(i)
   enddo

   read(in,*)nelem
   write(io,'(/,A,i5/)')' numero de elementos     :',nelem
   write(io,*)' conectividade dos elementos e propriedades'
   write(io,'(A,3(7x,A),5x,A)')'elemento','no1','no2','no3','material'
   allocate(kon(nne*nelem))
   allocate(imat(nelem))
   do j=1,nelem
      read(in,*)i,ic(1),ic(2),ic(3),imat(i)
      write (io,'(5i10)') i,ic(1),ic(2),ic(3),imat(i)
      n1=nne*(i-1)
      kon(n1+1)=ic(1)
      kon(n1+2)=ic(2)
      kon(n1+3)=ic(3)
   enddo

   read(in,*)nmat
   write(io,'(/,A,i5/)')' numero de materiais     :',nmat
   allocate(Ei(nmat))
   allocate(Poisson(nmat))
   do j=1,nmat
      read(in,*)i,Ei(i),Poisson(i)
      write(io,'(i10,2f20.5)')i,Ei(i),Poisson(i)
   enddo

   neq=nnode*ndf  !! calcula o numero total de incognitas
   allocate(Carga(neq))
   allocate(Desloc(neq))
   do i=1,neq
      Carga(i)=0.0 !! Zera vetor de cargas nodais
   enddo

   read(in,*)nln
   write(io,'(/,A,i5/)')' numero de nos carregados:',nln
   write(io,'(A,/7x,A,7x,A,8x,A)')'cargas nodais','no','px','py'
   do i=1,nln
      read(in,*) j,(w(k),k=1,ndf)
      write(io,'(i10,2f10.2)') j,(w(k),k=1,ndf)
      do k=1,ndf
         k1=ndf*(j-1)+k
         Carga(k1)=w(k)
      enddo
   enddo

   read(in,*)nbn
   write(io,'(//,A,i5//)')' numero de nos suportados:',nbn
   write(io,*)' dados das condicoes de contorno'
   write(io,'(23X,A,14X,A)')'status','valores prescritos'
   write(io,'(16X,A)')'(1:prescrito, 0:livre)'
   write(io,'(7x,A,10x,A,9x,A,16x,A,9x,A)')'no','u','v','u','v'
   allocate(istatus((ndf+1)*nbn))
   allocate(Prescrito(neq))
   do i=1,nbn
      read(in,*) j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      write(io,'(3i10,10x,2f10.4)') j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      L1=(ndf+1)*(i-1)+1
      L2=ndf*(j-1)
      istatus(L1)=j
      do k=1,ndf
         n1=L1+k
         n2=L2+k
         istatus(n1)=ic(k)
         Prescrito(n2)=w(k)
      enddo
   enddo

return
end subroutine Entrada_Dados



subroutine Imprime_Resultados()
use variaveis
implicit none
integer :: i,j,k1,k2

      write(io,'(//,A,/)')' RESULTADOS'

      write(io,'(A,/)')' deslocamentos nodais'
      write(io,'(8X,A,12X,A,14x,A)')'no','u','v'
      do i=1,nnode
        k1=ndf*(i-1)+1
        k2=k1+ndf-1
        write(io,'(i10,2E15.5)') i,(Desloc(j),j=k1,k2)
      enddo

      write(io,'(/,A,/)')' tensoes nos elementos'
      write(io,'(6X,A,4X,4(8X,A))')'NUMBER','R-STRESS ','Z-STRESS ','T-STRESS ','RZ-STRESS'
      do i=1,nelem
        write (io,'(6X,I4,5X,4e16.7)') i,stress(2,i),stress(1,i),stress(3,i),stress(4,i)
      enddo

   return
end subroutine Imprime_Resultados



subroutine Resultados_GiD()
use variaveis
implicit none
integer :: i,j,n1,k1,k2

   !open(UNIT=i33,file=arq1,status='UNKNOWN')
   !! Escreve arquivo da malha
   write(i33,*)'MESH "Malha" dimension 2 ElemType Triangle Nnode 3'
   write(i33,*)'Coordinates'
   do i=1,nnode
      write(i33,*)i,'  ',x(i),'  ',y(i)
   enddo
   write(i33,*)'end coordinates'
   write(i33,*)'Elements'
   do i=1,nelem
      n1=nne*(i-1)
      write(i33,'(5(2X,I8))') i,kon(n1+1),kon(n1+2),kon(n1+3),imat(i)
   enddo
   write(i33,*)'end elements'
   !close(UNIT=i33)

   !open(UNIT=i34,file=arq2,status='UNKNOWN')
   !! Escreve arquivo de resultados
   write(i34,*)'GiD Post Results File 1.0'
   write(i34,*)'GaussPoints "Malha_gauss" ElemType Triangle "Malha"'
   write(i34,*)'Number Of Gauss Points: 1'
   write(i34,*)'Nodes not included'
   write(i34,*)'Natural Coordinates: Internal'
   write(i34,*)'End gausspoints'
   write(i34,*)'Result "Desloc" "Load Analysis" ',1,' Vector OnNodes'
   write(i34,*)'ComponentNames "X-Desloc", "Y-Desloc"'
   write(i34,*)'Values'
   do i=1,nnode
      k1=ndf*(i-1)+1
      k2=k1+ndf-1
      write(i34,'(i10,2(2X,E15.5))')i,(Desloc(j),j=k1,k2)
   enddo
   write(i34,*)'End Values'
   write(i34,*)'Result "Gauss Stress" "Load Analysis" ',1,' ',&
              & 'PlainDeformationMatrix OnGaussPoints "Malha_gauss"'
   write(i34,*)'Values'
      do i=1,nelem
         write(i34,'(I3,1X,6(1X,E12.5))') i,stress(2,i),stress(1,i),stress(3,i),stress(4,i),0.0,0.0
      enddo
   write(i34,*)'End Values'
   !close(UNIT=i34)

   return
end subroutine Resultados_GiD



subroutine Fecha_Arquivos()
use variaveis

   close(UNIT=in)
   close(UNIT=io)
   close(UNIT=i33)
   close(UNIT=i34)

   return
end subroutine Fecha_Arquivos


end module Entrada_Saida



module Matriz_Elemento

contains

subroutine Monta_Matriz_Elemento(nel)
use variaveis
implicit none
integer :: nel,L,n1,n2,n3
integer :: i,j,k,IX,IZ,JX,JZ
real(4) :: D(4,4),XE(3,2),ZW(3),ZX(3),ZY(3)
real(4) :: C(6,6),DBA(4,6),A(6,6),B(4,6)
real(4) :: YM,PR,GG,THIRD,Z,ZK,DCON,DCON2,QLC,VOL

      YM=Ei(imat(nel))
      PR=Poisson(imat(nel))
      GG=YM/(2*(1.+PR))
      L=nne*(nel-1) !! nel = numero do elemento corrente
      n1=kon(L+1)   !! ponto nodal 1
      n2=kon(L+2)   !! ponto nodal 2
      n3=kon(L+3)   !! ponto nodal 1
      XE(1,1)=X(N1)
      XE(1,2)=Y(N1)
      XE(2,1)=X(N2)
      XE(2,2)=Y(N2)
      XE(3,1)=X(N3)
      XE(3,2)=Y(N3)
      do i=1,6
        do j=1,6
          RigElem(i,j)=0.0
        enddo
      enddo
      do J=1,6
        do I=1,4
          B(I,J)=0.0
          DB(I,J)=0.0
          DBA(I,J)=0.0
        enddo
        do I=1,6
          A(I,J)=0.0
          C(I,J)=0.0
        enddo
      enddo
      do J=1,4
        do I=1,4
          D(I,J)=0.0
        enddo
      enddo
      THIRD=1.0/3.0
      ORX=(XE(1,1)+XE(2,1)+XE(3,1))*THIRD
      ORY=(XE(1,2)+XE(2,2)+XE(3,2))*THIRD
      ZW(1)=XE(2,1)*XE(3,2)-XE(3,1)*XE(2,2)
      ZW(2)=XE(3,1)*XE(1,2)-XE(1,1)*XE(3,2)
      ZW(3)=XE(1,1)*XE(2,2)-XE(2,1)*XE(1,2)
      do I=1,3
        XE(I,1)=XE(I,1)-ORX
        XE(I,2)=XE(I,2)-ORY
      enddo
      ZX(1)=XE(2,2)-XE(3,2)
      ZX(2)=XE(3,2)-XE(1,2)
      ZX(3)=XE(1,2)-XE(2,2)
      ZY(1)=XE(3,1)-XE(2,1)
      ZY(2)=XE(1,1)-XE(3,1)
      ZY(3)=XE(2,1)-XE(1,1)
      ZK=XE(2,1)*XE(3,2)-XE(3,1)*XE(2,2)
      Z=3.0*ZK
!!    ELASTICITY MATRIX FOR A XI-SYMMETRIC CASE
      DCON=(YM*(1.-PR))/((1+PR)*(1.-2.*PR))
      DCON2=PR/(1.-PR)
      D(1,1)=DCON
      D(1,2)=DCON*DCON2
      D(1,3)=D(1,2)
      D(1,4)=0.0
      D(2,1)=D(1,2)
      D(2,2)=DCON
      D(2,3)=D(1,2)
      D(2,4)=0.0
      D(3,1)=D(1,3)
      D(3,2)=D(2,3)
      D(3,3)=DCON
      D(3,4)=0.0
      D(4,1)=D(1,4)
      D(4,2)=D(2,4)
      D(4,3)=D(3,4)
      D(4,4)=GG
!!    B MATRIX FOR AXI-SYMMETRIC CASE
      B(1,2)=ZY(1)
      B(1,4)=ZY(2)
      B(1,6)=ZY(3)
      B(2,1)=ZX(1)
      B(2,3)=ZX(2)
      B(2,5)=ZX(3)
      QLC=ORY/ORX
      B(3,1)=ZW(1)/ORX+ZX(1)+ZY(1)*QLC
      B(3,3)=ZW(2)/ORX+ZX(2)+ZY(2)*QLC
      B(3,5)=ZW(3)/ORX+ZX(3)+ZY(3)*QLC
      B(4,1)=ZY(1)
      B(4,2)=ZX(1)
      B(4,3)=ZY(2)
      B(4,4)=ZX(2)
      B(4,5)=ZY(3)
      B(4,6)=ZX(3)
      do I=1,4
        do J=1,6
          do K=1,4
            DB(I,J)=DB(I,J)+D(I,K)*B(K,J)/Z
          enddo
        enddo
      enddo
      VOL=3.14159265*Z*ORX
      if (vol.le.0.0) then
        write(io,*)'nel ',nel,'   volume ',vol
      endif
      do I=1,6
        do J=1,6
          do K=1,4
            RigElem(I,J)=RigElem(I,J)+B(K,I)*DB(K,J)*VOL/Z
          enddo
        enddo
      enddo

   return
end subroutine Monta_Matriz_Elemento



subroutine Calcula_Incognitas_Secundarias()
use variaveis
implicit none
integer :: i,j,k,L,jj,nel,NODLN
real(4) :: C(6,6),DBA(4,6)

      allocate(stress(4,nelem))
      do i=1,4
        do j=1,nelem
          stress(i,j)=0.0
        enddo
      enddo

      NODLN=1 ! Numero de casos de carregamento, sempre 1
      do NEL=1,nelem
        call Monta_Matriz_Elemento(nel)
        do J=1,NODLN
          do I=1,3
            L=NNE*(NEL-1)+1
            JJ=KON(L+I-1)
            C(2*I-1,J)=Desloc(2*JJ-1)
            C(2*I,J)=Desloc(2*JJ)
          enddo
        enddo
        do J=1,NODLN
          do I=1,4
            DBA(I,J)=0.0
            do K=1,6
              DBA(I,J)=DBA(I,J)+DB(I,K)*C(K,J)
            enddo
          enddo
        enddo
        do i=1,4
          stress(i,NEL)=DBA(i,1)
        enddo
      enddo

   return
end subroutine Calcula_Incognitas_Secundarias


end module Matriz_Elemento



module Matriz_Global

contains

subroutine Monta_Matriz_Global()
use variaveis
use Matriz_Elemento
implicit none
integer :: i,j,nel

    allocate(RigGlobal(neq,neq))
    do i=1,neq
       do j=1,neq
          RigGlobal(i,j)=0.0  !! Zera matriz global
       enddo
    enddo

    do nel=1,nelem
       call Monta_Matriz_Elemento(nel)
       call Armazena_Elemento_na_Global(nel)
    enddo

   return
end subroutine Monta_Matriz_Global



subroutine Armazena_Elemento_na_Global(nel)
use variaveis
implicit none
integer :: nel
integer :: i,j,i1,j1,i2,j2,n1,n2,k,L,jr,kr,jc,kc

   do i=1,nne
      n1=kon(nne*(nel-1)+i)
      i1=ndf*(i-1)
      j1=ndf*(n1-1)
      do j=1,nne
         n2=kon(nne*(nel-1)+j)
         i2=ndf*(j-1)
         j2=ndf*(n2-1)
         do k=1,ndf
            jr=i1+k
            kr=j1+k
            do L=1,ndf
               jc=i2+L
               kc=j2+L
               RigGlobal(kr,kc)=RigGlobal(kr,kc)+RigElem(jr,jc)
            enddo
         enddo
      enddo
   enddo

   return
end subroutine Armazena_Elemento_na_Global



subroutine Impoe_Condicoes_de_Contorno()
use variaveis
implicit none
integer :: L,no,k1,i,kr,j

   do L=1,nbn
      no=istatus((ndf+1)*(L-1)+1)
      k1=ndf*(no-1)
      do i=1,ndf
         if (istatus((ndf+1)*(L-1)+1+i).eq.1) then
            kr=k1+i
            do j=1,neq
               Carga(j)=Carga(j)-RigGlobal(kr,j)*Prescrito(kr)
               RigGlobal(kr,j)=0.0
               RigGlobal(j,kr)=0.0
            enddo
            RigGlobal(kr,kr)=1.0
            Carga(kr)=Prescrito(kr)
         endif
      enddo
   enddo

   return
end subroutine Impoe_Condicoes_de_Contorno



subroutine Sistema_Linear_Gauss(n,A,B,X,io)
implicit none
integer :: n,io
real(4) :: A(n,n),B(n),X(n)
integer :: i,j,k
real(4) :: c

   do k=1,(n-1)
      c=A(k,k)
      if (abs(c).lt.1.0E-7) then
         write(io,*)'Singularidade na linha ',k
         return
      endif
      do j=1,n
         A(k,j)=A(k,j)/c
      enddo
      B(k)=B(k)/c
      do i=(k+1),n
         c=A(i,k)
         do j=1,n
            A(i,j)=A(i,j)-c*A(k,j)
         enddo
         B(i)=B(i)-c*B(k)
      enddo
   enddo
   if (abs(A(n,n)).lt.1.0E-7) then
      write(io,*)'Singularidade na linha ',k
      return
   endif
   B(n)=B(n)/A(n,n)
   do i=(n-1),1,-1
      do j=(i+1),n
         B(i)=B(i)-A(i,j)*B(j)
      enddo
   enddo
   do i=1,n
      X(i)=B(i)
   enddo

   return
end subroutine Sistema_Linear_Gauss


end module Matriz_Global



program cst_axis95
use variaveis
use Entrada_Saida
use Matriz_Elemento
use Matriz_Global

implicit none
integer :: i,k1,k2,j
   call Abre_Arquivos()
   call Entrada_Dados()
   call Monta_Matriz_Global()
   call Impoe_Condicoes_de_Contorno()
   call Sistema_Linear_Gauss(neq,RigGlobal,Carga,Desloc,io)
   call Calcula_Incognitas_Secundarias()
   call Imprime_Resultados()
   call Resultados_GiD()
   call Fecha_Arquivos()

end program cst_axis95

