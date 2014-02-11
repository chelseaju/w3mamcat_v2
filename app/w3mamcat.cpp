//
//  w3mamcat.cpp
//  w3mamcat
//
//  Created by Chelsea Ju on 12/15/13.
//  Copyright (c) 2013 UCLA Biocybernetics. All rights reserved.
//

#include <cppcms/application.h>
#include <cppcms/service.h>
#include <cppcms/http_response.h>
#include <cppcms/url_dispatcher.h>
#include <cppcms/url_mapper.h>
#include <cppcms/applications_pool.h>
#include <cppcms/session_interface.h>
#include <iostream>
#include <stdlib.h>
#include "master.h"
#include "../model/model.h"


class w3mamcat :
public cppcms::application{
    
    public:
        w3mamcat(cppcms::service &srv) :
        cppcms::application(srv)  {
            dispatcher().assign("",&w3mamcat::intro,this);  
            mapper().assign("");  

            dispatcher().assign("/input",&w3mamcat::input,this);  
            mapper().assign("input","/input");  

            dispatcher().assign("/result",&w3mamcat::result,this);  
            mapper().assign("result","/result");  

            dispatcher().assign("/contact",&w3mamcat::contact,this);  
            mapper().assign("contact","/contact");  
        
            mapper().root("/w3mamcat");              
        }
    
    virtual void intro();
    virtual void contact();
    virtual void input();    
    virtual void result();
    void redirect(std::string);
    
    private:
    void setModelType(Model,std::string);
    void setOutputType(Model,std::string);
    void setTimeUnit(Model,std::string);
    void setVolUnit(Model,std::string);
    void setMassUnit(Model,std::string);
    void setEndoMU(Model,std::string);
//    void setBWType(Model,std::string);
//    void setErrorType(Model, std::string);

    std::vector<std::string> parseTFunction(Model);
    std::vector<std::string> parseIDParams(Model);
    std::string parseDerivedParams(Model);
    std::vector<std::string> parseCorrMatrix(Model);
    std::vector<std::string> parseMassVolume(Model); 
    std::vector<std::string> parseParameterBounds(Model);
    std::string parseIndividualBounds(float, float,float,float, std::string);
    std::string parseIndividualBoundsPlusmn(float, float,float,float, std::string, float[]);
    
};

int main(int argc,char ** argv)  
{  
    try {  
        cppcms::service srv(argc,argv);  
        srv.applications_pool().mount(  
                                      cppcms::applications_factory<w3mamcat>()  
                                      );  
        srv.run();  
    }  
    catch(std::exception const &e) {  
        std::cerr << e.what() << std::endl;  
    }  
}

void ini(content::master &c)  
{  
    c.intro_tab = "";
    c.input_tab = "";
    c.result_tab = "";
    c.contact_tab = "";
    c.help_tab = "";
}

void w3mamcat::intro()  
{  
    content::master c;
    ini(c);
    c.intro_tab = "active";
    render("intro",c);  
}  

void w3mamcat::contact()
{
    content::master c;
    c.contact_tab = "active";
    render("contact", c);
}

void w3mamcat::input()
{
    content::master c;
    ini(c);
    c.input_tab = "active";
             
    if(request().request_method()=="POST") {        
        c.user_input.load(context());  
        
        // clear data
        if(c.user_input.clear_form.value() != 0){
            c.user_input.clear();
            session().clear();            
        }
        else{        
            if(c.user_input.validate()){
                
                // if everything pass the validation, store the user input in session                c.mamcatModel.Reset();
                session().clear();
                
                session()["modelType"] = c.user_input.modelType.selected_id();
                session()["outputType"] = c.user_input.outputType.selected_id();
                session().set("nCompartment",  c.user_input.nCompartment.value());
                session().set("dose", c.user_input.dose.value());
                session()["mass"] = c.user_input.mass.selected_id();
                session()["time"] = c.user_input.time.selected_id();
                session()["volume"] = c.user_input.volume.selected_id();
                session()["endoMass"] = c.user_input.endoMass.selected_id();
                session()["sigDig"] = c.user_input.sigDig.selected_id();
                                
                
                if(c.user_input.bodyWeight.set()){
                    session().set("bodyWt", c.user_input.bodyWeight.value());
                }
                
                if(c.user_input.steadyStateValue.set()){
                    session().set("ssvalue", c.user_input.steadyStateValue.value());
                }
                
                for (int i = 1; i <= c.user_input.nCompartment.value(); i ++ )
                {
                    std::stringstream index;
                    index << i;
                    session().set("yt_Ainput_"+index.str(), c.user_input.yt_Ainput[i-1]->value());
                    session().set("yt_Linput_"+index.str(), c.user_input.yt_Linput[i-1]->value());                    
                }
               redirect("/result");
            }
            
            else
            {
                // make sure that the yt are not blocked
                int n = 10;
                if( c.user_input.nCompartment.valid())
                    n = c.user_input.nCompartment.value();

                for (int i = 1; i <= 10; i ++ )
                {
                    std::stringstream index;
                    index << i;
                    std::string id_A = "yt_Ainput_" + index.str();
                    std::string id_L = "yt_Linput_" + index.str();
                    
                    if(i <= n)
                    {
                        c.user_input.yt_Ainput[i-1] -> disabled(false);
                        c.user_input.yt_Linput[i-1] -> disabled(false);
                    }
                    else
                    {
                        c.user_input.yt_Ainput[i-1] -> disabled(true);
                        c.user_input.yt_Linput[i-1] -> disabled(true);
                    }
                    
                }
            }
                        
        }
    }    
    else{  // load the previously submitted information from session
        if(session().is_set("nCompartment")){
            int i = atoi(session()["nCompartment"].c_str());
            if(i > 0){
                c.user_input.nCompartment.value(i);
            }
        }
        if(session().is_set("modelType")){
            c.user_input.modelType.selected_id(session()["modelType"]);
        }
        if(session().is_set("outputType")){
            c.user_input.outputType.selected_id(session()["outputType"]);
        }
        if(session().is_set("dose")){
            double d = atof(session()["dose"].c_str());            
            if(d > 0){
                c.user_input.dose.value(d);
            }
        }
        
        if(session().is_set("mass")){
            c.user_input.mass.selected_id(session()["mass"]);
        }
        if(session().is_set("time")){
            c.user_input.time.selected_id(session()["time"]);
        }
        if(session().is_set("mass")){
            c.user_input.volume.selected_id(session()["volume"]);
        }
        if(session().is_set("endoMass")){
            c.user_input.endoMass.selected_id(session()["endoMass"]);
        }
        
        if(session().is_set("ssvalue")){
            double d = atof(session()["ssvalue"].c_str());
            if(d > 0){
                c.user_input.steadyStateValue.value(d);
            }
        }        
        if(session().is_set("bodywt")){
            double d = atof(session()["bodywt"].c_str());
            if(d > 0){
                c.user_input.bodyWeight.value(d);
            }
        }
        if(session().is_set("sigDig")){
            c.user_input.sigDig.selected_id(session()["sigDig"]);
        }

        for (int i = 1; i <= 10; i ++ )
        {
            std::stringstream index;
            index << i;
            std::string id_A = "yt_Ainput_" + index.str();
            std::string id_L = "yt_Linput_" + index.str();
            
            if(session().is_set(id_A))
                c.user_input.yt_Ainput[i-1] -> value(atof(session()[id_A].c_str()));
            if(session().is_set(id_L))
                c.user_input.yt_Linput[i-1] -> value(atof(session()[id_L].c_str()));
            
            if(i <= atoi(session()["nCompartment"].c_str()))
            {
                c.user_input.yt_Ainput[i-1] -> disabled(false);
                c.user_input.yt_Linput[i-1] -> disabled(false);
            }
            else
            {
                c.user_input.yt_Ainput[i-1] -> disabled(true);
                c.user_input.yt_Linput[i-1] -> disabled(true);
            }
            
        }
                
    }
    
    render("input",c);  

}

void w3mamcat::result()
{
    content::result c;
    ini(c);
    c.result_tab = "active";
    c.setmodel = false;
    
    if(session().is_set("modelType") && session().is_set("outputType") && session().is_set("time") &&
       session().is_set("mass") && session().is_set("volume") && session().is_set("endoMass") &&
       session().is_set("nCompartment") && session().is_set("dose"))
    {
        setModelType(c.mamcatModel, session()["modelType"]);
        setOutputType(c.mamcatModel, session()["outputType"]);
        setTimeUnit(c.mamcatModel, session()["time"]);
        setMassUnit(c.mamcatModel, session()["mass"]);
        setVolUnit(c.mamcatModel, session()["volume"]);
        setEndoMU(c.mamcatModel, session()["endoMass"]);
        c.mamcatModel.PoolNum() = atoi(session()["nCompartment"].c_str());
        c.mamcatModel.Dose() = atof(session()["dose"].c_str());

        // if compartment number if 1, Mam == Cat, set Model to be cat
        if(c.mamcatModel.PoolNum() == 1){
            c.mamcatModel.MType() = Mam;
        } 

        c.setmodel = true;
    }
       
    if(session().is_set("bodyWt")){
        c.mamcatModel.BodyWeight() = atof(session()["bodyWt"].c_str());
    }
    
    if(session().is_set("ssvalue")){
        c.mamcatModel.SSConc1() = atof(session()["ssvalue"].c_str());
    }

    for (int i = 1; i <= c.mamcatModel.PoolNum(); i ++ )
    {
        std::stringstream index;
        index << i;
        std::string id_A = "yt_Ainput_" + index.str();
        std::string id_L = "yt_Linput_" + index.str();
        
        if(session().is_set(id_A))
            c.mamcatModel.A(i) = atof(session()[id_A].c_str());
        if(session().is_set(id_L))
            c.mamcatModel.L(i) = atof(session()[id_L].c_str());        
    }
     
    // calculation
    if(c.setmodel == true){
        c.mamcatModel.Calculate();        
        c.tfFunction = parseTFunction(c.mamcatModel);
        c.idParameters = parseIDParams(c.mamcatModel);
        c.derivedParams = parseDerivedParams(c.mamcatModel);
        c.corrMatrix = parseCorrMatrix(c.mamcatModel);
        c.massVolume = parseMassVolume(c.mamcatModel);
        c.paramBound = parseParameterBounds(c.mamcatModel);
    
    }
    render("result",c);
            
}

void w3mamcat::redirect(std::string path)
{
    response().set_redirect_header(request().script_name() + path);
}



// private function
void w3mamcat::setModelType(Model m, std::string mtype)
{
    if(mtype == "mam")
        m.MType() = Mam;        
    else if (mtype == "cat")
        m.MType() = Cat;
}

void w3mamcat::setOutputType(Model m, std::string output)
{
    
    if(output == "con")
        m.OutputUnit() = Conc;
    else if (output == "mass")
        m.OutputUnit() = Mass;
}

void w3mamcat::setTimeUnit(Model m, std::string timeU)
{
    if(timeU == "sec")
        m.TimeUnit() = T_SEC;
    else if(timeU == "min")
        m.TimeUnit() = T_MIN;
    else if (timeU == "hr")
        m.TimeUnit() = T_HOUR;
    else if (timeU == "day")
        m.TimeUnit() = T_DAY;
}

void w3mamcat::setMassUnit(Model m, std::string bwt)
{
    if(bwt == "dose")
        m.MassUnit() = M_PDOSE;
    else if (bwt == "g")
        m.MassUnit() = M_G;
    else if (bwt == "kg")
        m.MassUnit() = M_KG;
    else if (bwt == "pd")
        m.MassUnit()= M_LB;
}


void w3mamcat::setVolUnit(Model m, std::string vol)
{    
    if(vol == "ml")
        m.VolUnit() = V_ML;
    else if (vol == "ul")
        m.VolUnit() = V_UL;
    else if (vol == "l")
        m.VolUnit() = V_L;
}

void w3mamcat::setEndoMU(Model m, std::string endo)
{    
    if(endo == "g")
        m.EndoMU() = ENDO_M_G;
    else if (endo == "kg")
        m.EndoMU() = ENDO_M_KG;
    else if (endo == "pd")
        m.EndoMU() = ENDO_M_LB;
}

// this code is copy from Soloman Russel's javascript PrintTFParams
std::vector<std::string> w3mamcat::parseTFunction(Model m)
{
    std::vector<std::string> v;
    std::string topString;
    std::string botString;
    std::string digit;
    
    v.resize(2);
    digit = "%." + session()["sigDig"]+ "f";
    
    topString = ""; // in the format of s^2 + beta(1)s + beta(2)
    botString = ""; // in the format of s^3 + alpha(1)s^2 + alpha(2)s + alpha(3)


    if(m.PoolNum() > 0){
        
        if(m.PoolNum() > 1){
            std::stringstream index;
            index << m.PoolNum();
            botString += "S<sup>" + index.str() + "</sup> + ";
        }
        else
            botString += "S + ";

        for (int i= m.PoolNum(); i >= 1; i --)
        {
            char alpha[10];
            char beta[10];
            sprintf(alpha, digit.c_str(), m.Alpha(i));
            sprintf(beta, digit.c_str(), m.Beta(i));

            topString += std::string(beta);
            botString += std::string(alpha);

            switch( i-1 )
            {
                case 0:
                    break;
                case 1:
                    topString += " S";
                    topString += " + ";
                    botString += " S";
                    botString += " + ";
                    break;
                default:
                    std::stringstream index;
                    index << i-1;
                    topString += " S<sup>"+index.str()+"</sup> + ";
                    botString += " S<sup>"+index.str()+"</sup> + ";
                    break;
            }
        }

    }
    
    v[0] = topString;
    v[1] = botString;
    return v;
}

// this code is copy from Soloman Russel's javascript PrintTFParams
std::vector<std::string> w3mamcat::parseIDParams(Model m)
{
    std::vector<std::string> v;
    std::string digit;
    digit = "%." + session()["sigDig"]+ "f";

    for( int i = 1; i <= m.PoolNum(); i++ ) {
        
        std::stringstream index, next_index;
        index << i;
        next_index << i+1;

        char k[10];
        sprintf(k, digit.c_str(), m.K(i,i));
        
        std::string kappa = "k<sub>" + index.str() + index.str() + "</sub> = " + k + " " + session()["time"] + "<sup>-1</sup>";
        
        v.push_back(kappa);

        
        if( i != m.PoolNum())
        {
            char gamma[10];
            sprintf(gamma, digit.c_str(), m.Gamma(i+1));

            std::string tmp;
            if(session()["modelType"] == "mam") // for mammilary model
                tmp = "k<sub>1" + next_index.str() + "</sub>k<sub>" + next_index.str() + "1</sub> = " + gamma + " " + session()["time"] + "<sup>-2</sup>"; 
            else
                tmp = "k<sub>" + index.str() + next_index.str() + "</sub>k<sub>"+ next_index.str() + index.str() +"</sub> = " + gamma + " " + session()["time"] + "<sup>-2</sup>"; 
            
            v.push_back(tmp);
        }                

    }
    return v;
}

// this code is copy from Soloman Russel's javascript PrintTFParams
std::string w3mamcat::parseDerivedParams(Model m)
{
    std::string message;
    std::string digit;
    digit = "%." + session()["sigDig"]+ "f";

    if(m.SSConc1() > 0){

        char pcr[10];
        char pr[10];
        char mrtLower[10];
        char mrtUpper[10];
        char vdLower[10];
        char vdUpper[10];

        sprintf(pcr, digit.c_str(), m.PCR());
        sprintf(pr, digit.c_str(), m.PR());
        sprintf(mrtLower, digit.c_str(), m.MRT(Lower));
        sprintf(mrtUpper, digit.c_str(), m.MRT(Upper));
        sprintf(vdLower, digit.c_str(), m.VD(Lower));
        sprintf(vdUpper, digit.c_str(), m.VD(Upper));

        
        message = "Pool 1 (Plasma?) Clearance Rate: PCR = " + std::string(pcr) + "(" + session()["volume"] + "/" + session()["time"] + ") <br/>\n"; 

        message += "Pool 1 Production Rate: PR =  PCR x C<sub>1</sub> = " + std::string(pr) + "(" + session()["endoMass"] + "/" + session()["time"] + ") <br/>\n"; 
        
        
        message += "Mean Residence Time in Whole System : " + std::string(mrtLower) + " < MRT (" + session()["time"] + ") < " + std::string(mrtUpper) + "<br/>\n";

        message += "Pool 1 Equivalent Distribution Volume : " + std::string(vdLower) + " < VD (" + session()["volume"] + ") < " + std::string(vdUpper) + "<br/>\n";
        
    }
    else
        message = "Not Avaliable. Steady State hasn't Defined yet.";
    
    return message;
}

// parsing the correlation matrix and store each element row-wise in a vector
// need further implementation
std::vector<std::string> w3mamcat::parseCorrMatrix(Model m)
{
    int size  = m.Corr_XDim() * m.Corr_YDim();
    std::vector<std::string> v;
    bool invalidValue = false;
    
    v.resize(size);
    
    if( size > 0 && m.SampleN() > 0){
    }
    
    return v;
    
}


// parsing the compartment masses and volumes in a vector
std::vector<std::string> w3mamcat::parseMassVolume(Model m)
{
    std::vector<std::string> v;
    
    for (int i=0; i < m.PoolNum(); i++) 
    {
        std::stringstream index;
        index << i+1;
        
        // volume
        float lowerUnV = m.V(i+1, Lower); // unconstrained
        float upperUnV = m.V(i+1, Upper);
        float lowerConV = m.ConV(i+1, Lower);  // constrained
        float upperConV = m.ConV(i+1, Upper);
                
        std::string parm = "V<sub>" + index.str() + "</sub> (" + session()["volume"] + ")";
        v.push_back(parseIndividualBounds(upperUnV, lowerUnV, upperConV, lowerConV, parm));

        // mass
        float lowerUnQ = m.Q(i+1, Lower); // unconstrained
        float upperUnQ = m.Q(i+1, Upper);
        float lowerConQ = m.ConQ(i+1, Lower);  // constrained
        float upperConQ = m.ConQ(i+1, Upper);

        parm = "Q<sub>" + index.str() + "</sub> (" + session()["mass"] + ")";        
        v.push_back(parseIndividualBounds(upperUnQ, lowerUnQ, upperConQ, lowerConQ, parm));

    }    
    
    // mass flux kijQi
    for(int i = 0; i < 3*m.PoolNum() - 2; i++ )
    {

        int one, two;
        std::string one_s, two_s;
        m.ItoIJ(i, &one, &two);
        
        float lowerUnFlux = m.Mflux(one, two, Lower); //unconstrained
        float upperUnFlux = m.Mflux(one, two, Upper);
        float lowerUnCon = m.ConMflux(one, two, Lower); // constrained
        float upperUnCon = m.ConMflux(one, two, Upper);

        // adjust the kij index for catenary model
        // replace k13 to k23; k31 to k32; k14 to k34; k41 to k43 ..etc
        if(session()["modelType"] == "cat"){
            
            if(abs(two - one) > 1)
            {
                if(one > two)
                    two = one - 1;
                else
                    one = two - 1;
            }                
        }

        one_s = std::string(1, one + '0');
        two_s = std::string(1, two + '0');
        std::string parm = "k<sub>" + one_s + two_s + "</sub>Q<sub>" + two_s + "</sub> (" + session()["endoMass"] + "/" + session()["time"] +  ")";        
        v.push_back(parseIndividualBounds(upperUnFlux, lowerUnFlux, upperUnCon, lowerUnCon, parm));
    }
    
    
    return v;
}

std::vector<std::string> w3mamcat::parseParameterBounds(Model m)
{
    std::vector<std::string> v;
    
    float MIN_DIFF = .000001;	// used for deal with problems of real number truncation
  
    v.resize(m.PoolNum());

    for (int i=0; i < m.PoolNum(); i++) 
    {
        int one, two;
        m.ItoIJ(i, &one, &two);
        std::string one_s = std::string(1, one + '0');
        std::string two_s = std::string(1, two + '0');
        
        // parameter bounds
        float lowerUn = m.K(one, two, Lower) ; // unconstrained
        float upperUn = m.K(one, two, Upper);
        float lowerCon = m.ConK(one, two, Lower);  // constrained
        float upperCon = m.ConK(one, two, Upper);
        
        float plusmn[4]; //[ unLCV, conLCV, conUCV, unUCV ];
        
        // unLCV, unUVC, conLCV, conUCV
        if(m.SampleN() > 0)
        {
            if( abs(lowerCon) < MIN_DIFF || abs(upperCon) < MIN_DIFF || abs(lowerUn) < MIN_DIFF || abs(upperUn) < MIN_DIFF ) 
            {
                if( abs(lowerUn) < MIN_DIFF)
                    plusmn[0] = 0;
                if( abs(lowerCon) < MIN_DIFF)
                    plusmn[1] = 0;
                if( abs(upperCon) < MIN_DIFF)
                    plusmn[2] = 0;
                if( abs(upperUn) < MIN_DIFF)
                    plusmn[3] = 0;
            }
            
            /* 
             *  TO BE IMPLEMENTED 
             *
             */
        }
        else
        {
            plusmn[0] = 0;
            plusmn[1] = 0;
            plusmn[2] = 0;
            plusmn[3] = 0;                   
        }
        
        std::string parm = "k<sub>" + one_s + two_s + "</sub> (" + session()["time"] + "<sup>-1</sup>)";
        v[i] =parseIndividualBoundsPlusmn(upperUn, lowerUn, upperCon, lowerCon, parm, plusmn);
    }
    
    
    return v;
}



// private function to be called by parseMassVolume
std::string w3mamcat::parseIndividualBounds(float upperUn, float lowerUn, float upperCon, float lowerCon, std::string parm)
{
    float MIN_DIFF = .000001;	// used for deal with problems of real number truncation

    std::string message;
    std::string red = "<tr style= 'color:red'>\n";
    std::string normal = "<tr>\n";
    std::string signs[4];
    char upUn[10];
    char loUn[10];
    char upCon[10];
    char loCon[10];

    std::string digit;
    digit = "%." + session()["sigDig"]+ "f";

    // convert float to string
    sprintf(upUn, digit.c_str(), upperUn);
    sprintf(loUn, digit.c_str(), lowerUn);
    sprintf(upCon, digit.c_str(), upperCon);
    sprintf(loCon, digit.c_str(), lowerCon);
    
    
    if(abs(lowerUn - lowerCon) < MIN_DIFF)
        signs[0] = "=";
    else
        signs[0] = "<";
    
    if( abs(lowerCon - upperCon) < MIN_DIFF )
    {
        signs[1] = "=";
        signs[2] = "=";
    }
    else
    {
        signs[1] = "<="; 
        signs[2] = "<=";
    }
    
    if( abs(upperCon - upperUn) < MIN_DIFF )
    {
        signs[3] = "=";
        message = normal;
    }
    else
    {
        signs[3] = "<";
        message = red;
    }
    
    
    message += "<td style='text-align:center'>" + std::string(loUn) + "</td>\n" +
                "<td style='text-align:center'>" + signs[0] + "</td>\n" +
                "<td style='text-align:center'>" + std::string(loCon) + "</td>\n" +
                "<td style='text-align:center'>" + signs[1] + "</td>\n" +
                "<td style='text-align:center'>" + parm + "</td>\n"+
                "<td style='text-align:center'>" + signs[2] + "</td>\n" +
                "<td style='text-align:center'>" + std::string(upCon) + "</td>\n" +
                "<td style='text-align:center'>" + signs[3] + "</td>\n" +
                "<td style='text-align:center'>" + std::string(upUn) + "</td>\n" +
                "</tr>\n";
    
    return message;
        
}

std::string w3mamcat::parseIndividualBoundsPlusmn(float upperUn, float lowerUn, float upperCon, float lowerCon, std::string parm, float plusmn[])
{
    float MIN_DIFF = .000001;	// used for deal with problems of real number truncation
    
    std::string message;
    std::string red = "<tr style= 'color:red'>\n";
    std::string normal = "<tr>\n";
    std::string signs[4];
    char upUn[10];
    char loUn[10];
    char upCon[10];
    char loCon[10];
    
    char unLCV[10];
    char conLCV[10];
    char unUCV[10];
    char conUCV[10];
    
    std::string digit;
    digit = "%." + session()["sigDig"]+ "f";
    
    // convert float to string
    sprintf(upUn, digit.c_str(), upperUn);
    sprintf(loUn, digit.c_str(), lowerUn);
    sprintf(upCon, digit.c_str(), upperCon);
    sprintf(loCon, digit.c_str(), lowerCon);

    sprintf(unLCV, digit.c_str(), plusmn[0]);
    sprintf(conLCV, digit.c_str(),plusmn[1]);
    sprintf(unUCV, digit.c_str(), plusmn[2]);
    sprintf(conUCV, digit.c_str(),plusmn[3]);
    
    
    if(abs(lowerUn - lowerCon) < MIN_DIFF)
        signs[0] = "=";
    else
        signs[0] = "<";
    
    if( abs(lowerCon - upperCon) < MIN_DIFF )
    {
        signs[1] = "=";
        signs[2] = "=";
    }
    else
    {
        signs[1] = "<="; 
        signs[2] = "<=";
    }
    
    if( abs(upperCon - upperUn) < MIN_DIFF )
    {
        signs[3] = "=";
        message = normal;
    }
    else
    {
        signs[3] = "<";
        message = red;
    }
    
    
    message += "<td style='text-align:center'>" + std::string(loUn) + "&plusmn;" + std::string(unLCV)+ "</td>\n" +
    "<td style='text-align:center'>" + signs[0] + "</td>\n" +
    "<td style='text-align:center'>" + std::string(loCon) + "&plusmn;" + std::string(conLCV) + "</td>\n" +
    "<td style='text-align:center'>" + signs[1] + "</td>\n" +
    "<td style='text-align:center'>" + parm + "</td>\n"+
    "<td style='text-align:center'>" + signs[2] + "</td>\n" +
    "<td style='text-align:center'>" + std::string(upCon) + "&plusmn;" + std::string(conUCV) + "</td>\n" +
    "<td style='text-align:center'>" + signs[3] + "</td>\n" +
    "<td style='text-align:center'>" + std::string(upUn) + "&plusmn;" + std::string(unUCV) + "</td>\n" +
    "</tr>\n";
    
    return message;
}



