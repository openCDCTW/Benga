import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';

export default class Dendrogram_view extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_dengrogram = this.query_dengrogram.bind(this);
	};

	query_dengrogram(){
		
		if(this.state.png_file == undefined){
			fetch('api/dendrogram/dendrogram/' + window.batchid, { method: 'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
				png_file: result.png_file, 
				pdf_file: result.pdf_file,
				svg_file: result.svg_file, 
				emf_file: result.emf_file, 
				newick_file: result.newick_file })));
			
		}else {
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.query_dengrogram();
		this.interval = setInterval(this.query_dengrogram, 10000);
	}



    render() {

    	if(this.state.png_file == undefined){

    		return(
    			<div>
    				<Paper square>
                        <Tabs value={false} centered>
                            <Tab label=" " disabled/>
                        </Tabs>
                    </Paper>
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font> Please hold on ... </font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <img src="https://svgshare.com/i/9N5.svg" />
                    </div>
                    <br />
                    <br />
                    <br />
                </div>
    		);
    	
    	}else{
    		return (
    				<div id="url">
    					<Paper square>
                        	<Tabs value={false} centered>
                            	<Tab label=" " disabled/>
                        	</Tabs>
                    	</Paper>
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <img src={this.state.svg_file} />
                        </div>
    					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <font>Download Dendrogram</font>
                        </div>

                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
        					
                            <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Png 
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
        					<a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Pdf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Svg
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                emf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                newick
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
    					<br />
    					<br />
    					<Link to="/upload_profile" style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">
                                <ReplyIcon />
                                &nbsp;&nbsp;
                                Back
                            </Button>
                        </Link>
                        <br />
    					<br />
                        <br />
                        <br />
    				</div>
        	);
    	}
        
    }
}
