{  
    "w3mamcat" : {
    	"media" : "/media",
    	"root" : "",
    },
    
    "service" : {  
        "api" : "http",  
        "port" : 8080  
    },  
    
    "http" : {  
        "script" :  "/w3mamcat" ,
        "rewrite" : [
                     { "regex" : "/media(/.*)?", "pattern" : "$0" },
                     ]
    },
    
    "file_server" :{
        "enable" : true,
        "document_root" : "."
    },
    
    "views" : {
        "paths" : [ "." ],  
        "skins" : [ "w3mamcat_skin" ]              
    },
    
    "session" : {  
        "expire" : "renew",  
        "timeout" : 604800,  
        "location" : "client",  
        "client" :      {  
            "hmac" :        "sha1",  
            "hmac_key" :    "3891bbf7f845fd4277008a63d72640fc13bb9a31"  
        }      
    }, 
    "security" : {  
        "csrf" : {  
            "enable" : true  
        }  
    },  

    
}  
