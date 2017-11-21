
  MODULE EnergyAbsorbed_Reflected


  USE Global
  USE EnclosureGeometry
  USE EnergyBundleLocation
  
  
  IMPLICIT NONE
  

  CONTAINS

   SUBROUTINE AbsorptionReflection

!******************************************************************************
!
! PURPOSE:			Checking whether the energy bundle absorbed or reflected 
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER		:: I,J, K, IOS, count
	REAL(Prec2)  :: R_absorbed
	   
!   R_absorbed		Random number generated is used to verify whether the 
!                   intercepted energy is absorbed or reflected by comparing
!                   it with surface absorptance

!    JDS 11-10-2006 added all the " if (WriteLogFile) then" blocks to control whether
!    or not a log file is written.  Also changed format statements to remove commas
!    so that RTVT could actually read the file.
!
    R_absorbed = Rand(6)     
!	   
IF (SpIndex .eq. 1) then    !Specular energy
     IF(R_absorbed < SpecReflec(SInter))Then
        NAEnergyS(SIndexR,SInter) = NAEnergyS(SIndexR,SInter) + 1   !RS: Total number of energy bundles absorbed
        TSpecA(SIndexR) = TSpecA(SIndexR) + 1    !RS: Total Number of energy bundles absorbed by surface
        
!   count the number of energy bundles absorbed and emitted
        IF (NCount .NE. 1) THEN !RS: Non-reflected rays
           NAEnergyWR(SIndexR,SInter)=NAEnergyWR(SindexR,Sinter)+1  !RS: Total Number of energy bundles absorbed without reflection
            
		  if (WriteLogFile) then
            Write(4,101, ADVANCE='YES')'P',SIndex,XLS(SIndex),YLS(SIndex),ZLS(SIndex),SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
          end if

101      Format(A1,2(' ',I2,3(' ',f6.3)),' ' '0')

        ELSE    !RS: Rays reflected or rereflected
         if (WriteLogFile) then
            Write(4,111,ADVANCE='YES')SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
         end if

111      Format(1(' ',I2,3(' ',f6.3)),' ' '0')
         
         NCount=0   !RS: Debugging: Resetting the counter
         
         TSpecR(SIndexR)=TSpecR(SIndexR)+1  !RS: Total number of reflected or reflected rays absorbed by surface

         NAEnergyR(SIndexR,SInter)=NAEnergyR(SIndexR,SInter)+1  !RS: Total number of reflected or rereflected rays absorbed
		 
        END IF		

         NTrials = NTrials + 1  !RS: Overall number of rays absorbed from the emitting surface
		   
        IF (SIndex ==1 .and. SInter == 2)Then
           count = count + 1	 

        ENDIf         

	    SIndex = SIndexR 
	    REF_IND = 0
	    Reflected = .False.   !RS: Setting the reflection flag to false since the array has now been absorbed
        TCountSpecR = 0 !RS:Resetting Reflection Flag
		          		     		
     ELSE   
        IF(SIndex == SIndexR .and. REF_IND == 0)Then    !RS: Reflected Rays
            
         NCount=1   !RS: Marking a reflection
                  
         TCountSpecR=1  !RS: Setting a flag

        if (WriteLogFile) then
           Write(4,112, ADVANCE='NO')'P',SIndex,XLS(SIndex),YLS(SIndex),ZLS(SIndex),SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
        end if
            
112     Format(A1,2(' ',I2,3(' ',f6.3)))

        ELSEIF(REF_IND == 1)Then    !RS: Rereflected rays
         
         TCountSpecR=2  !RS: Debugging: Setting a flag
         
           if (WriteLogFile) then
            Write(4,102,ADVANCE='NO')SInter, XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
           end if

102     format(' ',1(I2,3(' ',f6.3)))
        
        END IF
        
        OldSurface=SIndex  !RS: Keeping track of the emitting surface for each bounce
        SIndex=SInter
		 Reflected= .True.  !RS: Setting a flag for reflection
         REF_IND=1  !RS: A flag to say the ray has already been reflected at least once
         
     ENDIF	
Else    !Diffuse Energy
    IF(R_absorbed < Emit(SInter))Then
        NAEnergy(SIndexR,SInter) = NAEnergy(SIndexR,SInter) + 1
   	    TCOUNTA(SIndexR) = TCOUNTA(SIndexR) + 1

!   count the number of energy bundles absorbed and emitted

        IF (NCountd .NE. 1) THEN    !RS: Marking as a non-reflected ray
		  if (WriteLogFile) then
            Write(4,101, ADVANCE='YES')'P',SIndex,XLS(SIndex),YLS(SIndex),ZLS(SIndex),SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
          end if
        ELSE
		  if (WriteLogFile) then
            Write(4,111,ADVANCE='YES')SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
		  end if

		 NCountd=0  !RS: Resetting the reflection counter
         
        END IF		

         NTrialsd = NTrialsd + 1  
        
        IF (SIndex ==1 .and. SInter == 2)Then
           count = count + 1	 
        ENDIf         

	    SIndex = SIndexR 
	    REF_IND = 0
	    Reflected = .False.  !RS: Setting the reflection flag to false since the array has now been absorbed
		          		     		
     ELSE   
        IF(SIndex == SIndexR .and. REF_IND == 0)Then
         TCOUNTR(SIndexR) = TCOUNTR(SIndexR) + 1  

         NCountd = 1    !RS: Counting as a reflected surface
         
            if (WriteLogFile) then
                Write(4,112, ADVANCE='NO')'P',SIndex,XLS(SIndex),YLS(SIndex),ZLS(SIndex),SInter,XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
            end if
        ELSEIF(REF_IND == 1)Then
	     TCOUNTRR(SIndexR)= TCOUNTRR(SIndexR) + 1
         
           if (WriteLogFile) then
          Write(4,102,ADVANCE='NO')SInter, XP(SIndex,SInter),YP(SIndex,SInter),ZP(SIndex,SInter)
			end if

        END IF
        
        SIndex=SInter
		 Reflected= .True.
         REF_IND=1 
         
     ENDIF	
Endif

  END SUBROUTINE AbsorptionReflection
  END MODULE EnergyAbsorbed_Reflected




			 	
