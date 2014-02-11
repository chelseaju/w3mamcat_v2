function switch_on_off(itself)
{
	if(itself == true){
		e.disabled=false;
	}
	else{
		e.disabled=true;
	}
}

function switch_text(text,e)
{
	e.textContent = text;
}

function add_yt(e)
{
    for (var i = 1; i <= 10; i++)
    {
        var textfield_name_A = "ytAinput_" + i;
        var textfield_name_L = "ytLinput_" + i;
        
        var elem_A = document.getElementById(textfield_name_A);
        var elem_L = document.getElementById(textfield_name_L);
        
        if(i <= e.value)
        {
            elem_A.disabled = false;
            elem_L.disabled = false;
        }
        else
        {
            elem_A.value = '';
            elem_L.value = '';
            elem_A.disabled = true;
            elem_L.disabled = true;            
        }
    }
}

function add_error()
{
    var e = document.getElementsByName('errorType');
    var var_model = document.getElementsByName('var_model');
    var const_var = document.getElementById('constant_error');
    var var_disabled = false;
    var const_disabled = true;
    
    if(e[4].checked)
    {
        var_disabled= false;
        const_disabled = true;
    }
    else
    {
        var_disabled= true;
        const_disabled = false;
    }
    
    for (var i = 0; i < var_model.length; i++)
    {
        var_model[i].disabled = var_disabled;        
    }
    const_var.disabled = const_disabled;
    
}

function add_sample(e)
{
	var error = document.getElementById('error');
	var error_html = "<tr>\n" +
    "<td><b>Sample</b></td>\n" + 
    "<td><b>Time (sec)</b></td>\n"+
    "</tr>\n";
	for (var i=0; i< e.value; i++){
		error_html += "<tr>\n" +
        "<td><b>"+(i+1)+"</b></td>\n" + 
        "<td><input type='text' name='error' id='error_'"+(i+1)+" value=''></td>\n"+
        "</tr>\n";
	}
	error.innerHTML = error_html;
}

