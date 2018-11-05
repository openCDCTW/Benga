import React from 'react';
import ReactDOM from 'react-dom';

export default class Options extends React.Component {


    render(){
        return (
        	<form name="scheme">
        		<div>
                <label>Scheme: &nbsp;&nbsp;</label>
                    <select name="occurrence">
                    <option value="95">95</option>
                    <option value="90">90</option>
                    <option value="85">85</option>
                    <option value="80">80</option>
                    <option value="50">50</option>
                    </select>
                </div>
                <br />
                <div>
                <label>Database:&nbsp;&nbsp;</label>
                    <select name="database">
                    <option value="Salmonella_enterica">Salmonella_enterica</option>
                    <option value="Vibrio_cholerae">Vibrio cholerae</option>
                    <option value="Listeria_monocytogenes">Listeria monocytogenes</option>
                    <option value="Campylobacter_jejuni">Campylobacter jejuni</option>
                    </select>
                </div>
                <br />
                <div>
                <label>Email:&nbsp;&nbsp;&nbsp; </label>
                <input name="email" type="text" placeholder="Input your Email here"/>
                </div>
        	</form>
        );
    }
}