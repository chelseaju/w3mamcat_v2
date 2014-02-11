//
//  master.h
//  w3mamcat
//
//  Created by Chelsea Ju on 12/15/13.
//  Copyright (c) 2013 UCLA Biocybernetics. All rights reserved.
//


#include <cppcms/view.h>
#include <cppcms/form.h>
#include <iostream>
#include "../model/model.h"

namespace content  {
    
    struct input_form: public cppcms::form{

        // model type
        cppcms::widgets::radio modelType;
        cppcms::widgets::numeric<int> nCompartment;
        
        // output type
        cppcms::widgets::radio outputType;
        cppcms::widgets::numeric<double> steadyStateValue;
        
        // unit
        cppcms::widgets::select mass;  
        cppcms::widgets::select time;  
        cppcms::widgets::select volume;  
        cppcms::widgets::select endoMass;
        
        // other info
        cppcms::widgets::select sigDig;
        cppcms::widgets::numeric<double> dose;
        cppcms::widgets::numeric<double> bodyWeight;
        
        // error model
        cppcms::widgets::numeric<int> sampleSize;
        cppcms::widgets::radio errorType;
        cppcms::widgets::numeric<double> con_error;
        cppcms::widgets::numeric<double> var_b;
        cppcms::widgets::numeric<double> var_c;
        cppcms::widgets::numeric<double> var_d;
        
        
        // yt input
        cppcms::form yt_Ainput_form;
        cppcms::form yt_Linput_form;
        std::vector< cppcms::widgets::numeric<double>* > yt_Ainput; 
        std::vector< cppcms::widgets::numeric<double>* > yt_Linput; 
                
        // submit
        cppcms::widgets::submit submit;

        // clear data
        cppcms::widgets::submit clear_form;
                
        input_form(){
            nCompartment.message("Number of Compartments");
            steadyStateValue.message("Steady State Value (optional)");
            
            mass.message("Mass");
            time.message("Time");
            volume.message("Volume");
            endoMass.message("Endo Mass");
            
            sigDig.message("Sig. Digit for Display");
            dose.message("Dose");
            bodyWeight.message("Body Wt Normalization (optional)");
            
            sampleSize.message("Sample Size");
            con_error.message("Constant (CV or SD or VAR) = ");
            var_b.message("B = ");
            var_c.message("C = ");
            var_d.message("D = ");

            errorType.name("errorType");
            con_error.id("constant_error");
            var_b.name("var_model");
            var_c.name("var_model");
            var_d.name("var_model");
            
            
            submit.value("Run Analysis");
            clear_form.value("Clear Data");
            
            // adding field to the form
            add(modelType);
            add(nCompartment);
            
            add(outputType);
            add(steadyStateValue);
            
            add(mass);
            add(time);
            add(volume);
            add(endoMass);
            
            add(sigDig);
            add(dose);
            add(bodyWeight);
            
            add(sampleSize);
            add(errorType);
            add(con_error);
            add(var_b);
            add(var_c);
            add(var_d);
            
            add(submit);
            add(clear_form);
                        
            add(yt_Linput_form);
            add(yt_Ainput_form);
            for (int i = 0; i < 10; i++)
            {
                std::stringstream index;
                index << i+1;
                cppcms::widgets::numeric<double> * a = new cppcms::widgets::numeric<double>();
                a -> message(index.str());
                a -> id("ytAinput_" + index.str());
                a -> disabled(true);
                yt_Ainput_form.attach(a);
                yt_Ainput.push_back(a);

                cppcms::widgets::numeric<double> * l = new cppcms::widgets::numeric<double>();
                l -> id("ytLinput_" + index.str());
                l -> disabled(true);
                yt_Linput_form.attach(l);
                yt_Linput.push_back(l);
            
            }
            
            
            // adding selection values
            modelType.add("Mammillary", "mam");
            modelType.add("Catenary", "cat");
            
            outputType.add("Concentration (gm/ml)", "con");
            outputType.add("Mass (gm)", "mass");
            
            mass.add("dose", "dose");
            mass.add("gram", "g");
            mass.add("kilogram", "kg");
            mass.add("pound", "pd");
            
            time.add("seconds", "sec");
            time.add("minutes", "min");
            time.add("hours", "hr");
            time.add("days", "day");
            
            volume.add("milliliter", "ml");
            volume.add("microliter", "ul");
            volume.add("liter", "l");
            
            endoMass.add("grams", "g");
            endoMass.add("kilograms", "kg");
            endoMass.add("pound", "pd");
            
            for (int i = 1; i < 10; i++){                
                std::string t = std::string(1, i+ '0');
                sigDig.add(t, t);
            }
            
            errorType.add("Constant CV (%)", "cv");
            errorType.add("Constant SD", "sd");
            errorType.add("Constant VAR", "var");
            errorType.add("Variable CV (%) to Specify", "var_cv");
            errorType.add("VAR = B + Cz^D", "var_eqn");
            
            con_error.disabled(true);
            var_b.disabled(true);
            var_c.disabled(true);
            var_d.disabled(true);
            
            // adding error message
            modelType.error_message("Must select model type");
            nCompartment.error_message("Must specify a valid number (1-10)");
            outputType.error_message("Must select output type");
            dose.error_message("Must specify dose");
                        
            // adding constraints
            modelType.non_empty();
            nCompartment.non_empty();
            nCompartment.low(1);
            nCompartment.high(10);
            outputType.non_empty();
            dose.non_empty();
        }
        
        
        virtual bool validate()
        {
            if(!form::validate()) 
                return false;
            
            // checking the input constraints
            if(nCompartment.value() > 0 && nCompartment.value() <= 10){
                int n = nCompartment.value();
                for (int i = 0; i < 10; i ++)
                {
                    if( i < n && (!yt_Ainput[i]->set() || yt_Ainput[i]->value() <= 0 ))
                    {
                        yt_Ainput[i] -> valid(false);
                        yt_Ainput[i] -> disabled(true);
                        yt_Ainput[i] -> error_message("Must be positive");
                        return false;                    
                    }
                    
                    if( i < n && (!yt_Linput[i]->set() || yt_Linput[i]->value() >= 0 ))
                    {
                        yt_Linput[i] -> valid(false);
                        yt_Linput[i] -> disabled(true);
                        yt_Linput[i] -> error_message("Must be negative");
                        return false;                    
                    }
                    
                    if(i >= n)
                    {
                        yt_Linput[i] -> clear();
                        yt_Ainput[i] -> clear();
                        yt_Linput[i] -> disabled(true);
                        yt_Ainput[i] -> disabled(true);
                    }
                }
            }
            
            
            if(dose.value() <= 0){
                dose.valid(false);
                dose.error_message("Must be Positive");
                return false;
            }
            
            if(steadyStateValue.set() && steadyStateValue.value() <= 0){
                steadyStateValue.valid(false);
                steadyStateValue.error_message("Must be positive");
                return false;
            }
            
            if(bodyWeight.set() && bodyWeight.value() <= 0){
                bodyWeight.valid(false);
                bodyWeight.error_message("Must be positive");
                return false;
            }
            
            if(mass.set() && mass.selected_id() == "dose" && dose.set() && dose.value() != 100)
            {
                dose.error_message("Dose needs to be 100 since A(i) units are %Dose");
                dose.valid(false);
                return false;
            }
            
            if(errorType.set() && errorType.selected_id() == "var_eqn" && (!var_b.set() || !var_c.set() || !var_d.set()))
            {
                if(!var_b.set())
                {
                    var_b.error_message("Missing data");
                    var_b.valid(false);                    
                }
                if(!var_c.set())
                {
                    var_c.error_message("Missing data");
                    var_c.valid(false);                    
                }
                if(!var_d.set())
                {
                    var_d.error_message("Missing data");
                    var_d.valid(false);                    
                }
                
                return false;
            }
            
            if(errorType.set() && errorType.selected_id() != "var_eqn" && !con_error.set())
            {
                con_error.error_message("Missing data");
                con_error.valid(false);
                return false;
            }
            
            
           return true;
        }
        
    };
    
    
    struct master : public cppcms::base_content{
        std::string intro_tab;
        std::string input_tab;
        std::string result_tab;
        std::string contact_tab;
        std::string help_tab;
        Model mamcatModel;
        input_form user_input;
    };
    
    struct result : public master {        
        bool setmodel; // a flag to indicate whether a mamcat model is created        
        std::vector<std::string> tfFunction; // a vector to store transfer function
                                               // it should only contain two elements, the first one is the nominators (alphas) and the seoncd one is the denominators (betas)
        std::vector<std::string> idParameters; // contains identifiable parameters + their values
        std::string derivedParams ; // whole organism derived parameters
        std::vector<std::string> corrMatrix; 
        std::vector<std::string> massVolume; // compartment masses and volumes
        std::vector<std::string> paramBound; // parameter constrains
    };
    
    
    
    
}