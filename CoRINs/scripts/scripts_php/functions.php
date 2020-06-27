<?php
   
    
    if(isset($_SESSION['count_name_session'])){
               
        $count_name_session = $_SESSION['count_name_session'];
        $count_name_session++; 
        $count_name_session = $_SESSION['count_name_session'] = $count_name_session;            
        
             
        mkdir('../../result_web/'.$count_name_session.'/', 0777, true);           
        
        //echo '<br /><a href="../../site/index.html">Pagina FELIPE PROJETO</a>';
    } else {
        $count_name_session = $_SESSION['count_name_session'] = 0;              
    }
        
?>