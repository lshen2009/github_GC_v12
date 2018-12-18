Code.GC12: the currently working directory
Code.GC12_Nov4: Delete the last species and see if I can make the code run. But it does not work.
Code.GC12_Nov7: Add LSHENA in the Standard.eqn and then delete it using KPP. Now the code works,
                but not in an elegant way.
Code.GC12_Nov8: Now the code works well based on Code.GC12_Nov7, but in an elegant way.
Code.GC12_Nov9: The number of Functions that needs editting is reduced to 5. I test after deleting 6 species.
Code.GC12_Nov10: I test two mechanisms, one for whole chem in the troposphere and the other for nonVoCs in the Statosphere. Now the simulation time can reduce 20%.
Code.GC12_Nov12: I defined var_selected and var_deleted in the Global variable. Now the simulation time can reduce 20%.
Code.GC12_Nov13: Now the gckpp_Integrator can adapt to varying dimensions of fast and slow species.
Code.GC12_Nov15: Now the code can deal with both fast and slow species. But the simulation speed doesnot improve.
Code.GC12_Nov17: Based on Code.GC12_Nov15, but the code becomes more elegant.
Code.GC12_Nov18: Based on Code.GC12_Nov17, I have added all 10 types of chemistry regimes. But the speed does not improve.
Code.GC12_Nov20: The code works for 10 tyeps of chemistry regimes.
Code.GC12_Nov21: The code works for 15 tyeps of chemistry regimes.
Code.GC12_Nov21_afternoon: The code works for 15 tyeps of chemistry regimes. I have reduced the number of input parameters.
Code.GC12_Nov25: The code works for 15 tyeps of chemistry regimes. I consider both fast and slow species in FUN and Jac_SP functions.
Code.GC12_Dec04: The code works for 15 tyeps of chemistry regimes. I consider both fast and slow species in FUN and Jac_SP functions. But the model may be wrong in clean regions and in the stratosphere.
Code.GC12_Dec16: The code works for 15 tyeps of chemistry regimes. I consider both fast and slow species in FUN and Jac_SP functions. But the model may be wrong in clean regions and in the stratosphere. This is the version after my Dec05 group meeting.I find the Het Function may cause discrepancies. 
Code.GC12_Dec18: Now I am trying to find the day-night line and assign this as the full-chem regime.
