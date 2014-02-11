#line 1 "view/master.tmpl"
#include "app/master.h" 
#line 3 "view/master.tmpl"
namespace w3mamcat_skin {
	#line 4 "view/master.tmpl"
	struct master :public cppcms::base_view
	#line 4 "view/master.tmpl"
	{
	#line 4 "view/master.tmpl"
		content::master &content;
	#line 4 "view/master.tmpl"
		master(std::ostream &_s,content::master &_content): cppcms::base_view(_s),content(_content)
	#line 4 "view/master.tmpl"
		{
	#line 4 "view/master.tmpl"
		}
		#line 6 "view/master.tmpl"
		virtual void title() {
			#line 7 "view/master.tmpl"
			out()<<"\n"
				"    ";
			#line 7 "view/master.tmpl"
			out()<<cppcms::locale::translate("W3MAMCAT - Mamcat Parameter Estimator");
			#line 8 "view/master.tmpl"
			out()<<"\n"
				"";
		#line 8 "view/master.tmpl"
		} // end of template title
		#line 10 "view/master.tmpl"
		virtual void header() {
			#line 11 "view/master.tmpl"
			out()<<"\n"
				"    <title>";
			#line 11 "view/master.tmpl"
			title();
			#line 15 "view/master.tmpl"
			out()<<"</title>\n"
				"    <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n"
				"    <link href=\"/media/style.css\" rel=\"stylesheet\" type=\"text/css\" />\n"
				"    <script src=\"/media/javascript.js\"></script>\n"
				"";
		#line 15 "view/master.tmpl"
		} // end of template header
		#line 17 "view/master.tmpl"
		virtual void page_content() {
			#line 19 "view/master.tmpl"
			out()<<"\n"
				"    Override  Me\n"
				"";
		#line 19 "view/master.tmpl"
		} // end of template page_content
		#line 21 "view/master.tmpl"
		virtual void render() {
			#line 26 "view/master.tmpl"
			out()<<"\n"
				"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\"\n"
				"    \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n"
				"<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
				"<head>\n"
				"    ";
			#line 26 "view/master.tmpl"
			header();
			#line 43 "view/master.tmpl"
			out()<<"\n"
				"</head>\n"
				"\n"
				"<body>\n"
				"\n"
				"<div id=\"main\" class=\"box\">\n"
				"\n"
				"  <!--- LOGO --->\n"
				"  <div id=\"header\">\n"
				"    <h1 id=\"logo\"><strong>W<sup>3</sup></strong>MAMCAT v2.0</h1>\n"
				"    <hr class=\"noscreen\" />\n"
				"  </div>\n"
				"\n"
				"  <!--- NAVIGATION TAB --->\n"
				"  <div id=\"tabs\" class=\"noprint\">\n"
				"    <h3 class=\"noscreen\">Navigation</h3>\n"
				"    <ul class=\"box\">\n"
				"      <li id='";
			#line 43 "view/master.tmpl"
			out()<<cppcms::filters::escape(content.intro_tab);
			#line 43 "view/master.tmpl"
			out()<<"'><a href='";
			#line 43 "view/master.tmpl"
			content.app().mapper().map(out(),"");
			#line 44 "view/master.tmpl"
			out()<<"'>Home</a></li>\n"
				"      <li id='";
			#line 44 "view/master.tmpl"
			out()<<cppcms::filters::escape(content.input_tab);
			#line 44 "view/master.tmpl"
			out()<<"'><a href='";
			#line 44 "view/master.tmpl"
			content.app().mapper().map(out(),"input");
			#line 45 "view/master.tmpl"
			out()<<"'>Inputs</a></li>\n"
				"      <li id='";
			#line 45 "view/master.tmpl"
			out()<<cppcms::filters::escape(content.result_tab);
			#line 45 "view/master.tmpl"
			out()<<"'><a href='";
			#line 45 "view/master.tmpl"
			content.app().mapper().map(out(),"result");
			#line 46 "view/master.tmpl"
			out()<<"'>Results</a></li>\n"
				"      <li id='";
			#line 46 "view/master.tmpl"
			out()<<cppcms::filters::escape(content.help_tab);
			#line 48 "view/master.tmpl"
			out()<<"'><a href=\"\">Help</a></li>\n"
				"      <li><a href=\"http://www.biocyb.cs.ucla.edu/w3mamcat/expert/1.htm\" target=\"_blank\">MAMCAT v1.0</a></li>\n"
				"      <li id='";
			#line 48 "view/master.tmpl"
			out()<<cppcms::filters::escape(content.contact_tab);
			#line 48 "view/master.tmpl"
			out()<<"'><a href='";
			#line 48 "view/master.tmpl"
			content.app().mapper().map(out(),"contact");
			#line 54 "view/master.tmpl"
			out()<<"'>Contacts</a></li>\n"
				"    </ul>\n"
				"    <hr class=\"noscreen\" />\n"
				"  </div>\n"
				"\n"
				"  <div id=\"page\" class=\"box\">\n"
				"    ";
			#line 54 "view/master.tmpl"
			page_content();
			#line 68 "view/master.tmpl"
			out()<<"  \n"
				" </div>\n"
				" <hr class=\"noscreen\" />\n"
				" \n"
				" <div id=\"footer\">\n"
				"    <p id=\"copyright\">created by Chelsea Ju & Sepideh Mazrouee</br>\n"
				"    &copy; 2013 Prof. Joe DiStefano III, <a href=\"http://www.biocyb.cs.ucla.edu/\">UCLA Biocybernetics</a> | All Right Reserved</p>\n"
				" </div>\n"
				"\n"
				" </div>\n"
				"</body>\n"
				"\n"
				"</html>\n"
				"\n"
				"";
		#line 68 "view/master.tmpl"
		} // end of template render
	#line 69 "view/master.tmpl"
	}; // end of class master
#line 70 "view/master.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/contact.tmpl"
namespace w3mamcat_skin {
	#line 2 "view/contact.tmpl"
	struct contact :public master
	#line 2 "view/contact.tmpl"
	{
	#line 2 "view/contact.tmpl"
		content::master &content;
	#line 2 "view/contact.tmpl"
		contact(std::ostream &_s,content::master &_content): master(_s,_content),content(_content)
	#line 2 "view/contact.tmpl"
		{
	#line 2 "view/contact.tmpl"
		}
		#line 4 "view/contact.tmpl"
		virtual void page_content() {
			#line 40 "view/contact.tmpl"
			out()<<"  \n"
				"\n"
				"\t<div id=\"content\">\n"
				"        <div class=\"article\">\n"
				"            <h2>Contacts</h2>\n"
				"            <div class=\"contact\">\t\t\n"
				"                <b>Joseph J. DiStefano III</b><br/>\n"
				"                Professor of Computer Science, Medicine, and Biomedical Engineering<br/>\n"
				"                UCLA<br/>\n"
				"                4711 Boelter Hall<br/>\n"
				"                Los Angeles, CA 90095-1596 <br/>\n"
				"                (310) 825-7482 <br/>\n"
				"                <a href=\"mailto:joed@cs.ucla.edu\">joed@cs.ucla.edu</a><br/>\n"
				"            </div>\n"
				"\n"
				"            <div class=\"contact\">\n"
				"                <b>Chelsea Ju</b><br/>\n"
				"                Graduate Student <br/>\n"
				"                UCLA Computer Science Department<br/> \n"
				"                <a href=\"mailto:chelseaju@cs.ucla.edu\">chelseaju@cs.ucla.edu</a><br/>\n"
				"            </div>\n"
				"\n"
				"            <div class=\"contact\">\n"
				"                <b>Sepideh Mazrouee</b><br/>\n"
				"                Graduate Student <br/>\n"
				"                UCLA Computer Science Department<br/>\n"
				"                <a href=\"mailto:sepideh@cs.ucla.edu\">sepideh@cs.ucla.edu</a><br/>\n"
				"            </div>\n"
				"\n"
				"            <div class=\"contact\">\n"
				"                <b>UCLA Biocybernetics Laboratory </b><br/>\n"
				"                <a href=\"mailto:http://www.biocyb.cs.ucla.edu\">http://www.biocyb.cs.ucla.edu</a><br/>\n"
				"            </div>\t\n"
				"\n"
				"        </div>\n"
				"\t</div>\n"
				"";
		#line 40 "view/contact.tmpl"
		} // end of template page_content
	#line 41 "view/contact.tmpl"
	}; // end of class contact
#line 42 "view/contact.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/intro.tmpl"
namespace w3mamcat_skin {
	#line 2 "view/intro.tmpl"
	struct intro :public master
	#line 2 "view/intro.tmpl"
	{
	#line 2 "view/intro.tmpl"
		content::master &content;
	#line 2 "view/intro.tmpl"
		intro(std::ostream &_s,content::master &_content): master(_s,_content),content(_content)
	#line 2 "view/intro.tmpl"
		{
	#line 2 "view/intro.tmpl"
		}
		#line 4 "view/intro.tmpl"
		virtual void page_content() {
			#line 41 "view/intro.tmpl"
			out()<<"  \n"
				"\n"
				"\t<div id=\"content\">\n"
				"\t\t<div class=\"article\">\t\t\n"
				"\t\t\t<h2>Parameter Estimator</h2>\n"
				"\t\t\t<p><b>1. MODELING:</b> A linear n-compartment MAMMILLARY and/or CATENARY model can be formulated and fully quantified by tranformation from a known MULTIEXPONENTIAL MODEL (sum of n exponentials). Solutions are limited to models with both input and output in compartment 1, and n must be no greater than 10. <br/>\n"
				"\t\t\tNOTE: A 2-step modeling process is assumed. The multiexponential model is typically obtained from kinetic time-series data using a data fitting program like DIMSUM+, another expert system available from the Biocybernetics Laboratory.\n"
				"\t\t\t</p>\n"
				"\t\t\t<p><b>2. MODEL DISCRIMINATION/DISTINGUISHABILITY:</b> An expert subsystem is included for determining whether a mammillary can be distiguished form a catenary linear compartmental model, each with n compartments, but input and output may be in any (same of different) compartments. We assume both can be fitted to a single output data set; the questions is whether the two different structures can be distinguished via this data. The theory is fully developed pedagogically, along with specific solutions.\n"
				"\t\t\t</p> \n"
				"        </div>\n"
				"        <hr class=\"noscreen\" />\n"
				"\t\t<div class=\"article\">\n"
				"\t\t\t<h2>Previous Work</h2>\n"
				"\t\t\t<p>W<sup>3</sup>MAMCAT is the online version of MAMCAT+ Version 1.0, principally programmed by Hsiau-Te Su. The mathematical computations are done using the same code as the original program. The expert system for distinguishing mammillary and catenary models was originally programmed by Anh-Ngoc B. Kuhn De Chizelle. This version (v2.0) is an extention of Solomon Russell's work. The web application is designed to follow the MVC framework using CppCMS. The expert system is hosted on the server side, and user's inputs are stored in session, and parsed securely to the server. Results are then transfered back to the browser of client side.\n"
				"            </p>   \n"
				"            \n"
				"\t\t</div>\n"
				"        <hr class=\"noscreen\" />\n"
				"\t\t<div class=\"resource\">\n"
				"\t\t\t<h2>References</h2>\n"
				"\t\t\t<ul>\n"
				"\t\t\t\t<li>Box, Don. <i>Essential COM</i>. Addison-Wesley Professional, 1998.</li>\n"
				"\t\t\t\t<li>Chappell, David. <i>Understanding ActiveX and OLE: a guide for developers and managers</i>. Microsoft Press, 1996.</li>\n"
				"\t\t\t\t<li>Distefano III, J. J. \"Complete parameter bounds and quasiidentifiability conditions for a class of unidentifiable linear systems.\" <i>Mathematical Biosciences</i> 65.1 (1983): 51-68.</li>\n"
				"\t\t\t\t<li>E.M. Landaw, B. C. Chen, J. J. DiStefano, III. \"An Algorithm for the Identifiable Parameter Combinations of the General.\" <i>Mathematical Biosciences</i> 72.2 (1982):199-212.</li>\n"
				"\t\t\t\t<li>Chen, Benjamin Chao-Min, Elliot M. Landaw, and Joseph J. DiStefano. \"Algorithms for the identifiable parameter combinations and parameter bounds of unidentifiable catenary compartmental models.\" <i>Mathematical biosciences</i> 76.1 (1985): 59-68.</li>\n"
				"\t\t\t\t<li>Distefano, Joseph J., Benjamin C. Chen, and Elliot M. Landaw. \"Pool size and mass flux bounds and quasiidentifiability relations for catenary models.\" <i>Mathematical Biosciences</i> 88.1 (1988): 1-14.</li>\n"
				"\t\t\t\t<li>Lindell, Robert, Joseph J. Distefano, and Elliot M. Landaw. \"Statistical variability of parameter bounds for< i> n</i>-pool unidentifiable mammillary and catenary compartmental models.\" <i>Mathematical biosciences</i> 91.2 (1988): 175-199.</li>\n"
				"\t\t\t\t<li>Kuhn de Chizelle, Anh-Ngoc B., and Joseph J. DiStefano. \"MAMCAT: an expert system for distinguishing between mammillary and catenary compartmental models.\" <i>Computers in biology and medicine</i> 24.3 (1994): 189-204.</li>\n"
				"\t\t\t\t<li>Vicini, Paolo, Hsiao-Te Su, and Joseph J. Distefano Iii. \"Identifiability and interval identifiability of mammillary and catenary compartmental models with some known rate constants.\" <i>Mathematical biosciences</i> 167.2 (2000): 145-161.</li>\n"
				"                <li>Russell, Solomon, and Joseph J. DiStefano III. \"W< sup> 3</sup> MAMCAT: A world wide web based tool for mammillary and catenary compartmental modeling and expert system distinguishability.\" <i>Computer methods and programs in biomedicine</i> 83.1 (2006): 34-42.</li>\n"
				"                <li>\"CppCMS - The C Web Development Framework.\" CppCMS. Publisher, Feb. 2013, Web. http://cppcms.com/wikipp </li>\n"
				"\n"
				"\t\t\t</ul>\n"
				"\t\t</div>\n"
				"    </div>\n"
				"";
		#line 41 "view/intro.tmpl"
		} // end of template page_content
	#line 42 "view/intro.tmpl"
	}; // end of class intro
#line 43 "view/intro.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/input.tmpl"
namespace w3mamcat_skin {
	#line 2 "view/input.tmpl"
	struct input :public master
	#line 2 "view/input.tmpl"
	{
	#line 2 "view/input.tmpl"
		content::master &content;
	#line 2 "view/input.tmpl"
		input(std::ostream &_s,content::master &_content): master(_s,_content),content(_content)
	#line 2 "view/input.tmpl"
		{
	#line 2 "view/input.tmpl"
		}
		#line 4 "view/input.tmpl"
		virtual void page_content() {
			#line 9 "view/input.tmpl"
			out()<<"  \n"
				"\n"
				"\t<div id=\"content\">\n"
				"        <div class=\"io\">\n"
				"        \n"
				"        <form action=\"\" method=\"post\">\t";
			#line 9 "view/input.tmpl"
			out() << "<input type=\"hidden\" name=\"_csrf\" value=\"" << content.app().session().get_csrf_token() <<"\" >\n";
			#line 10 "view/input.tmpl"
			out()<<"\n"
				"            ";
			#line 10 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.clear_form).render(_form_context); }
			#line 11 "view/input.tmpl"
			out()<<"\n"
				"            ";
			#line 11 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.submit).render(_form_context); }
			#line 21 "view/input.tmpl"
			out()<<"\n"
				"\n"
				"            <table>\n"
				"                <tr>\n"
				"                <td>\n"
				"                <fieldset>\n"
				"                    <legend>Model Type</legend>\n"
				"                    <table>\n"
				"                    <tr>\n"
				"                        <td>\n"
				"                        ";
			#line 21 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.modelType).render(_form_context); }
			#line 30 "view/input.tmpl"
			out()<<"\n"
				"                        </td>\n"
				"                        <td>\n"
				"                            <img src='/media/mam_but.gif' style=\"height:10px\"><br/>\n"
				"                            <img src='/media/cat_but.gif' style=\"height:10px\">\n"
				"                        </td>\n"
				"                    </tr>\n"
				"                    <tr>\n"
				"                        <td colspan=2>\n"
				"                        ";
			#line 30 "view/input.tmpl"
			if(!(content.user_input.nCompartment.valid())) {
				#line 32 "view/input.tmpl"
				out()<<"\n"
					"                            <span class=\"cppcms_form_error\">\n"
					"                            ";
				#line 32 "view/input.tmpl"
				out()<<cppcms::filters::escape(content.user_input.nCompartment.error_message());
				#line 34 "view/input.tmpl"
				out()<<"\n"
					"                            </span>\n"
					"                        ";
			#line 34 "view/input.tmpl"
			} // endif
			#line 37 "view/input.tmpl"
			out()<<"\n"
				"                        \n"
				"                        <span class=\"cppcms_form_input\">\n"
				"                        ";
			#line 37 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 37 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::first_part);
			#line 37 "view/input.tmpl"
			    (content.user_input.nCompartment).render_input(_form_context);
			#line 37 "view/input.tmpl"
			}
			#line 39 "view/input.tmpl"
			out()<<"\n"
				"                        onChange = \"add_yt(this);\"\n"
				"                        ";
			#line 39 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 39 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::second_part);
			#line 39 "view/input.tmpl"
			    (content.user_input.nCompartment).render_input(_form_context);
			#line 39 "view/input.tmpl"
			}
			#line 52 "view/input.tmpl"
			out()<<"\n"
				"                        - Compartments\n"
				"                        </span>\n"
				"                        </td>\n"
				"                    </tr>\n"
				"                    </table>\n"
				"                </fieldset>\n"
				"                </td>\n"
				"\n"
				"                <td>\n"
				"                <fieldset>\n"
				"                    <legend>Units</legend>\n"
				"                    <table>\n"
				"                        ";
			#line 52 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.mass).render(_form_context); }
			#line 53 "view/input.tmpl"
			out()<<"\n"
				"                        ";
			#line 53 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.time).render(_form_context); }
			#line 54 "view/input.tmpl"
			out()<<" \n"
				"                        ";
			#line 54 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.volume).render(_form_context); }
			#line 55 "view/input.tmpl"
			out()<<"\n"
				"                        ";
			#line 55 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.endoMass).render(_form_context); }
			#line 71 "view/input.tmpl"
			out()<<"\n"
				"                    </table>\n"
				"                </fieldset>\n"
				"                </td>\n"
				"                \n"
				"                <td rowspan=\"2\">\n"
				"                <fieldset>\n"
				"                    <legend><i>y<sub>(t)</sub>= ∑ A<sub>i</sub>e<sup>L<sub>i</sub>t</sup></i></legend>\n"
				"                    <table>\n"
				"                    <tr style=\"height:20px\">\n"
				"                        <td>\n"
				"                        <table>\n"
				"                            <tr>\n"
				"                                <th align=\"center\">#</th>\n"
				"                                <th align=\"center\">A (% dose/ml)</th>\n"
				"                            </tr>\n"
				"                            ";
			#line 71 "view/input.tmpl"
			if((content.user_input.yt_Ainput_form).begin()!=(content.user_input.yt_Ainput_form).end()) {
				#line 72 "view/input.tmpl"
				out()<<"\n"
					"                                ";
				#line 72 "view/input.tmpl"
				for(CPPCMS_TYPEOF((content.user_input.yt_Ainput_form).begin()) a_ptr=(content.user_input.yt_Ainput_form).begin(),a_ptr_end=(content.user_input.yt_Ainput_form).end();a_ptr!=a_ptr_end;++a_ptr) {
				#line 72 "view/input.tmpl"
				CPPCMS_TYPEOF(*a_ptr) &a=*a_ptr;
					#line 73 "view/input.tmpl"
					out()<<"\n"
						"                                ";
					#line 73 "view/input.tmpl"
					{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (a).render(_form_context); }
					#line 74 "view/input.tmpl"
					out()<<"\n"
						"                                ";
				#line 74 "view/input.tmpl"
				} // end of item
				#line 75 "view/input.tmpl"
				out()<<"                                \n"
					"                            ";
			#line 75 "view/input.tmpl"
			}
			#line 85 "view/input.tmpl"
			out()<<"\n"
				"                        </table>\n"
				"                        </td>\n"
				"                        \n"
				"                        <td>\n"
				"                        <table>\n"
				"                        <tr>\n"
				"                            <th></th>\n"
				"                            <th align=\"center\">L (1/sec)</th>\n"
				"                        </tr>\n"
				"                            ";
			#line 85 "view/input.tmpl"
			if((content.user_input.yt_Linput_form).begin()!=(content.user_input.yt_Linput_form).end()) {
				#line 86 "view/input.tmpl"
				out()<<"\n"
					"                                ";
				#line 86 "view/input.tmpl"
				for(CPPCMS_TYPEOF((content.user_input.yt_Linput_form).begin()) l_ptr=(content.user_input.yt_Linput_form).begin(),l_ptr_end=(content.user_input.yt_Linput_form).end();l_ptr!=l_ptr_end;++l_ptr) {
				#line 86 "view/input.tmpl"
				CPPCMS_TYPEOF(*l_ptr) &l=*l_ptr;
					#line 86 "view/input.tmpl"
					{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (l).render(_form_context); }
				#line 86 "view/input.tmpl"
				} // end of item
				#line 87 "view/input.tmpl"
				out()<<"                                \n"
					"                            ";
			#line 87 "view/input.tmpl"
			}
			#line 105 "view/input.tmpl"
			out()<<"\n"
				"                        </table>\n"
				"                        </td>\n"
				"\n"
				"                    </tr>\n"
				"                    \n"
				"                    </table>\n"
				"\t\t\t\t</fieldset>\n"
				"                </td>                \n"
				"                </tr>\n"
				"            \n"
				"                <tr>\n"
				"                <td>\n"
				"                <fieldset>\n"
				"                    <legend>Output Type</legend>\n"
				"                    <table>\n"
				"                        <tr>\n"
				"                            <td colspan = 2>\n"
				"                                ";
			#line 105 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.outputType).render(_form_context); }
			#line 108 "view/input.tmpl"
			out()<<"\n"
				"                            </td>\n"
				"                        </tr>                            \n"
				"                        ";
			#line 108 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.steadyStateValue).render(_form_context); }
			#line 116 "view/input.tmpl"
			out()<<"                    \n"
				"                    </table>\n"
				"                </fieldset>\n"
				"                </td>\n"
				"                <td>\n"
				"                <fieldset>\n"
				"                    <legend>Other Info</legend>\n"
				"                    <table>\n"
				"                        ";
			#line 116 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.sigDig).render(_form_context); }
			#line 117 "view/input.tmpl"
			out()<<"\n"
				"                        ";
			#line 117 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.dose).render(_form_context); }
			#line 118 "view/input.tmpl"
			out()<<"\n"
				"                        ";
			#line 118 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_table); (content.user_input.bodyWeight).render(_form_context); }
			#line 133 "view/input.tmpl"
			out()<<"\n"
				"                    </table>\n"
				"                </fieldset>\n"
				"                </td>\n"
				"\n"
				"                </tr>\n"
				"\n"
				"                <tr>\n"
				"                    <td colspan = 2>\n"
				"                    <fieldset>\n"
				"                    <legend>Error Model (optional)</legend>\n"
				"                    <table>\n"
				"                    <tr>\n"
				"                        <td>\n"
				"                        Sample Size\n"
				"                        ";
			#line 133 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 133 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::first_part);
			#line 133 "view/input.tmpl"
			    (content.user_input.sampleSize).render_input(_form_context);
			#line 133 "view/input.tmpl"
			}
			#line 135 "view/input.tmpl"
			out()<<"\n"
				"                            onChange = \"add_sample(this)\"\n"
				"                        ";
			#line 135 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 135 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::second_part);
			#line 135 "view/input.tmpl"
			    (content.user_input.sampleSize).render_input(_form_context);
			#line 135 "view/input.tmpl"
			}
			#line 145 "view/input.tmpl"
			out()<<"\n"
				"                        </td>\n"
				"                        <td rowspan = 2>                                           \n"
				"                            <table id=\"error\">\n"
				"                            </table>\n"
				"                        </td>\n"
				"                    </tr>\n"
				"                    <tr>\n"
				"                        <td>\n"
				"                            Error e <sub>(tk)</sub> has: <br/>\n"
				"                            ";
			#line 145 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 145 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::first_part);
			#line 145 "view/input.tmpl"
			    (content.user_input.errorType).render_input(_form_context);
			#line 145 "view/input.tmpl"
			}
			#line 147 "view/input.tmpl"
			out()<<"\n"
				"                                onChange = \"add_error()\"\n"
				"                            ";
			#line 147 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html);
			#line 147 "view/input.tmpl"
			    _form_context.widget_part(cppcms::form_context::second_part);
			#line 147 "view/input.tmpl"
			    (content.user_input.errorType).render_input(_form_context);
			#line 147 "view/input.tmpl"
			}
			#line 153 "view/input.tmpl"
			out()<<"\n"
				"                         </td>\n"
				"                    </tr>\n"
				"                    \n"
				"                    <tr>\n"
				"                        <td>\n"
				"                            ";
			#line 153 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.con_error).render(_form_context); }
			#line 155 "view/input.tmpl"
			out()<<"<br/>\n"
				"                            VAR Model: <br/>\n"
				"                            ";
			#line 155 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.var_b).render(_form_context); }
			#line 156 "view/input.tmpl"
			out()<<" <br/>\n"
				"                            ";
			#line 156 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.var_c).render(_form_context); }
			#line 157 "view/input.tmpl"
			out()<<" <br/>\n"
				"                            ";
			#line 157 "view/input.tmpl"
			{ cppcms::form_context _form_context(out(),cppcms::form_flags::as_html,cppcms::form_flags::as_space); (content.user_input.var_d).render(_form_context); }
			#line 168 "view/input.tmpl"
			out()<<" <br/>\n"
				"                        </td>\n"
				"                    </tr>\n"
				"                    </table>\n"
				"                    </fieldset>                    \n"
				"                    </td>\n"
				"                </tr>\n"
				"            </table>\n"
				"        </form>  \n"
				"        </div>\n"
				"    </div>\n"
				"";
		#line 168 "view/input.tmpl"
		} // end of template page_content
	#line 169 "view/input.tmpl"
	}; // end of class input
#line 170 "view/input.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/result.tmpl"
namespace w3mamcat_skin {
	#line 2 "view/result.tmpl"
	struct result :public master
	#line 2 "view/result.tmpl"
	{
	#line 2 "view/result.tmpl"
		content::result &content;
	#line 2 "view/result.tmpl"
		result(std::ostream &_s,content::result &_content): master(_s,_content),content(_content)
	#line 2 "view/result.tmpl"
		{
	#line 2 "view/result.tmpl"
		}
		#line 3 "view/result.tmpl"
		virtual void page_content() {
			#line 12 "view/result.tmpl"
			out()<<"  \n"
				"    <div id=\"content\">\n"
				"\t\t<div class=\"io\">\n"
				"\t\t\t<table>\n"
				"\t\t\t\t<tr>\n"
				"\t\t\t\t\t<td>\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Transfer Function</legend>\n"
				"                        \n"
				"                        ";
			#line 12 "view/result.tmpl"
			if(content.setmodel) {
				#line 15 "view/result.tmpl"
				out()<<"\n"
					"                        <table style=\"padding:15px\";>\n"
					"                            <tr>\n"
					"                                <td style=\"border-bottom:thin solid \"><center><i>";
				#line 15 "view/result.tmpl"
				out()<<cppcms::filters::raw(content.tfFunction.at(0));
				#line 18 "view/result.tmpl"
				out()<<" </i></center></td>\n"
					"                            </tr>\n"
					"                            <tr>\n"
					"                                <td><center><i>";
				#line 18 "view/result.tmpl"
				out()<<cppcms::filters::raw(content.tfFunction.at(1));
				#line 21 "view/result.tmpl"
				out()<<"</i></center></td>\n"
					"                            </tr>\n"
					"                        </table>\n"
					"                        ";
			#line 21 "view/result.tmpl"
			} // endif
			#line 28 "view/result.tmpl"
			out()<<"\n"
				"                        \n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t\t<td>\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Identifiable Parameter Combinations</legend>\n"
				"                        ";
			#line 28 "view/result.tmpl"
			if(content.setmodel) {
				#line 30 "view/result.tmpl"
				out()<<"\n"
					"                        <table>\n"
					"                            ";
				#line 30 "view/result.tmpl"
				if((content.idParameters).begin()!=(content.idParameters).end()) {
					#line 31 "view/result.tmpl"
					out()<<"\n"
						"                                ";
					#line 31 "view/result.tmpl"
					for(CPPCMS_TYPEOF((content.idParameters).begin()) combo_ptr=(content.idParameters).begin(),combo_ptr_end=(content.idParameters).end();combo_ptr!=combo_ptr_end;++combo_ptr) {
					#line 31 "view/result.tmpl"
					CPPCMS_TYPEOF(*combo_ptr) &combo=*combo_ptr;
						#line 34 "view/result.tmpl"
						out()<<"\n"
							"                                <tr>\n"
							"                                    <td>\n"
							"                                        <center>";
						#line 34 "view/result.tmpl"
						out()<<cppcms::filters::raw(combo);
						#line 37 "view/result.tmpl"
						out()<<"</center>\n"
							"                                    </td>\n"
							"                                </tr>\n"
							"                                ";
					#line 37 "view/result.tmpl"
					} // end of item
					#line 38 "view/result.tmpl"
					out()<<"\n"
						"                            ";
				#line 38 "view/result.tmpl"
				}
				#line 40 "view/result.tmpl"
				out()<<"\n"
					"                        </table>\n"
					"                        ";
			#line 40 "view/result.tmpl"
			} // endif
			#line 48 "view/result.tmpl"
			out()<<"\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>  \n"
				"                </tr>\n"
				"\t\t\t\t<tr>\n"
				"\t\t\t\t\t<td colspan=\"2\">\t\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t<legend>Parameter Constraints</legend>\n"
				"                    ";
			#line 48 "view/result.tmpl"
			if(!(content.paramBound.empty())) {
				#line 62 "view/result.tmpl"
				out()<<"\n"
					"                        <table>\n"
					"                            <tr>\n"
					"                                <td><font color='#7f000b'><center>Unconstrained &plusmn; CV</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Constrained &plusmn; CV</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>P</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Constrained &plusmn; CV</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Unconstrained &plusmn; CV</center></font></td>\n"
					"                            </tr>\n"
					"                            \n"
					"                            ";
				#line 62 "view/result.tmpl"
				if((content.paramBound).begin()!=(content.paramBound).end()) {
					#line 63 "view/result.tmpl"
					out()<<"\n"
						"                                ";
					#line 63 "view/result.tmpl"
					for(CPPCMS_TYPEOF((content.paramBound).begin()) bound_ptr=(content.paramBound).begin(),bound_ptr_end=(content.paramBound).end();bound_ptr!=bound_ptr_end;++bound_ptr) {
					#line 63 "view/result.tmpl"
					CPPCMS_TYPEOF(*bound_ptr) &bound=*bound_ptr;
						#line 64 "view/result.tmpl"
						out()<<"\n"
							"                                ";
						#line 64 "view/result.tmpl"
						out()<<cppcms::filters::raw(bound);
						#line 65 "view/result.tmpl"
						out()<<"</center>\n"
							"                                ";
					#line 65 "view/result.tmpl"
					} // end of item
					#line 66 "view/result.tmpl"
					out()<<"\n"
						"                            ";
				#line 66 "view/result.tmpl"
				}
				#line 68 "view/result.tmpl"
				out()<<"\n"
					"                        </table>\n"
					"                    ";
			#line 68 "view/result.tmpl"
			} // endif
			#line 77 "view/result.tmpl"
			out()<<"\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t</tr>\n"
				"                \n"
				"                <tr>\n"
				"\t\t\t\t\t<td colspan = 2>                    \n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Compartment Masses and Volumes</legend>\n"
				"                        ";
			#line 77 "view/result.tmpl"
			if(!(content.massVolume.empty())) {
				#line 91 "view/result.tmpl"
				out()<<"\n"
					"                        <table>\n"
					"                            <tr>\n"
					"                                <td><font color='#7f000b'><center>Unconstrained</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Constrained</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>P</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Constrained</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>&lt =</center></font></td>\n"
					"                                <td><font color='#7f000b'><center>Unconstrained</center></font></td>\n"
					"                            </tr>\n"
					"                            \n"
					"                            ";
				#line 91 "view/result.tmpl"
				if((content.massVolume).begin()!=(content.massVolume).end()) {
					#line 92 "view/result.tmpl"
					out()<<"\n"
						"                                ";
					#line 92 "view/result.tmpl"
					for(CPPCMS_TYPEOF((content.massVolume).begin()) bound_ptr=(content.massVolume).begin(),bound_ptr_end=(content.massVolume).end();bound_ptr!=bound_ptr_end;++bound_ptr) {
					#line 92 "view/result.tmpl"
					CPPCMS_TYPEOF(*bound_ptr) &bound=*bound_ptr;
						#line 93 "view/result.tmpl"
						out()<<"\n"
							"                                ";
						#line 93 "view/result.tmpl"
						out()<<cppcms::filters::raw(bound);
						#line 94 "view/result.tmpl"
						out()<<"</center>\n"
							"                                ";
					#line 94 "view/result.tmpl"
					} // end of item
					#line 95 "view/result.tmpl"
					out()<<"\n"
						"                            ";
				#line 95 "view/result.tmpl"
				}
				#line 97 "view/result.tmpl"
				out()<<"\n"
					"                        </table>\n"
					"                        ";
			#line 97 "view/result.tmpl"
			} // endif
			#line 106 "view/result.tmpl"
			out()<<"\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t</tr>\n"
				"\n"
				"\t\t\t\t<tr>\n"
				"\t\t\t\t\t<td colspan=\"2\">\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Whole Organism (Derived) Parameters</legend>\n"
				"                        ";
			#line 106 "view/result.tmpl"
			if(content.setmodel) {
				#line 108 "view/result.tmpl"
				out()<<"\n"
					"                        <center>\n"
					"                            ";
				#line 108 "view/result.tmpl"
				out()<<cppcms::filters::raw(content.derivedParams);
				#line 110 "view/result.tmpl"
				out()<<"\n"
					"                        </center>\n"
					"                        ";
			#line 110 "view/result.tmpl"
			} // endif
			#line 118 "view/result.tmpl"
			out()<<"\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t</tr>\n"
				"\t\t\t\t<tr>\n"
				"\t\t\t\t\t<td colspan=\"2\">\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Correlation Matrix</legend>\n"
				"                        ";
			#line 118 "view/result.tmpl"
			if(content.setmodel) {
				#line 119 "view/result.tmpl"
				out()<<"\n"
					"                           ";
				#line 119 "view/result.tmpl"
				if(content.corrMatrix.empty()) {
					#line 121 "view/result.tmpl"
					out()<<"\n"
						"                           <center><font color=red> NA: not enough sample size </font></center>\n"
						"                           ";
				#line 121 "view/result.tmpl"
				} // endif
				#line 123 "view/result.tmpl"
				out()<<"\n"
					"                           \n"
					"                        ";
			#line 123 "view/result.tmpl"
			} // endif
			#line 140 "view/result.tmpl"
			out()<<"\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t</tr>\n"
				"\n"
				"\t\t\t\t<tr>\n"
				"\t\t\t\t\t<td colspan=\"2\">\n"
				"\t\t\t\t\t<fieldset>\n"
				"\t\t\t\t\t\t<legend>Model Display</legend>\n"
				"\t\t\t\t\t</fieldset>\n"
				"\t\t\t\t\t</td>\n"
				"\t\t\t\t</tr>\n"
				"\n"
				"\t\t\t</tbody></table>\t\t\t\t\t\n"
				"\t\t</div>\n"
				"\t</div>\n"
				"\n"
				"";
		#line 140 "view/result.tmpl"
		} // end of template page_content
	#line 141 "view/result.tmpl"
	}; // end of class result
#line 142 "view/result.tmpl"
} // end of namespace w3mamcat_skin
#line 3 "view/master.tmpl"
namespace w3mamcat_skin {
#line 70 "view/master.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/contact.tmpl"
namespace w3mamcat_skin {
#line 42 "view/contact.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/intro.tmpl"
namespace w3mamcat_skin {
#line 43 "view/intro.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/input.tmpl"
namespace w3mamcat_skin {
#line 170 "view/input.tmpl"
} // end of namespace w3mamcat_skin
#line 1 "view/result.tmpl"
namespace w3mamcat_skin {
#line 142 "view/result.tmpl"
} // end of namespace w3mamcat_skin
#line 143 "view/result.tmpl"
namespace {
#line 143 "view/result.tmpl"
 cppcms::views::generator my_generator; 
#line 143 "view/result.tmpl"
 struct loader { 
#line 143 "view/result.tmpl"
  loader() { 
#line 143 "view/result.tmpl"
   my_generator.name("w3mamcat_skin");
#line 143 "view/result.tmpl"
   my_generator.add_view<w3mamcat_skin::master,content::master>("master",true);
#line 143 "view/result.tmpl"
   my_generator.add_view<w3mamcat_skin::contact,content::master>("contact",true);
#line 143 "view/result.tmpl"
   my_generator.add_view<w3mamcat_skin::intro,content::master>("intro",true);
#line 143 "view/result.tmpl"
   my_generator.add_view<w3mamcat_skin::input,content::master>("input",true);
#line 143 "view/result.tmpl"
   my_generator.add_view<w3mamcat_skin::result,content::result>("result",true);
#line 143 "view/result.tmpl"
    cppcms::views::pool::instance().add(my_generator);
#line 143 "view/result.tmpl"
 }
#line 143 "view/result.tmpl"
 ~loader() {  cppcms::views::pool::instance().remove(my_generator); }
#line 143 "view/result.tmpl"
} a_loader;
#line 143 "view/result.tmpl"
} // anon 
