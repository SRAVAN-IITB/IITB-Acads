library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity FSMInternet is
	port(b:in std_logic_vector(3 downto 0); clock:in std_logic; 
		led: out std_logic_vector(4 downto 0));
end FSMInternet;

architecture bhv of FSMInternet is
	
	type state is (idle,stud,rese, prof, admin, trans); 
	signal count:integer:=0;
	signal present_state, next_state: state:=idle;
	
	
	
	begin
		only_process:process(clock, b, present_state)
		begin
			case present_state is
				when idle=> ---------------------------------------------------------idle state
					led<="00000";
					
					if(b="0000") then 
						next_state<= idle;
					
					elsif (b="0001") then
						next_state<=stud;
					
					elsif (b="0011" or b="0010") then
						next_state<=rese;
					
					elsif (b(3)='0' and b(2)='1') then
						next_state<=prof;
					
					else
						next_state<=admin;
					end if;
					
				when stud=> ---------------------------------------------------------student state
					led<="00001";
					
					if(b="0000") then 
						next_state<= idle;
					
					elsif (b="0001") then
						next_state<=stud;
					
					else
						next_state<=trans;
					end if;
					
				when rese=> ---------------------------------------------------------researcher state
					led<="00011";
					
					if(b="0000" or b="0001") then 
						next_state<= idle;
					
					elsif (b="0011" or b="0010") then
						next_state<=rese;
					
					else
						next_state<=trans;
					end if;
					
				when prof=> ---------------------------------------------------------Professor state
					led<="00111";
					
					if(b(3) ='1') then 
						next_state<=trans;
					
					elsif (b(3)='0' and b(2)='1') then
						next_state<=prof;
					
					else
						next_state<=idle;
					end if;
				
				when admin=> ---------------------------------------------------------Admin state
					led<="01111";
					
					if(b(3)='0') then 
						next_state<= idle;
					
					else
						next_state<=admin;
					end if;
					
				when trans=> ---------------------------------------------------------Transition state
					led<="11111";
					
					if ((b="0011" or b="0010")) then
						next_state<=rese;
						
					
					elsif (b(3)='0' and b(2)='1') then
						next_state<=prof;
						

					elsif (b(3)='1') then
						next_state<=admin;
						
					else
						next_state<=trans;
						
					end if;
			end case;

			if(clock='1' and clock' event) then
				present_state<=next_state;
				
			else
				present_state<=present_state;
			end if;
		
		end process;
		
end architecture bhv;