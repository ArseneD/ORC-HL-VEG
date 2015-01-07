! Obsolete parallel tools : Barrier

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_parallel/tools_para.f90 $ 
!< $Date: 2014-02-20 17:26:48 +0100 (Thu, 20 Feb 2014) $
!< $Author: josefine.ghattas $
!< $Revision: 1925 $
!-

MODULE tools_para
!-
  USE mod_orchidee_para_var, ONLY : MPI_COMM_ORCH
!-
#include "src_parallel.h"
!-
CONTAINS

  SUBROUTINE barrier_para()
#ifdef CPP_PARA
    CALL MPI_BARRIER(MPI_COMM_ORCH,ierr)
#endif
  END SUBROUTINE barrier_para

END MODULE tools_para
