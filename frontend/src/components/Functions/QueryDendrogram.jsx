import React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';

const styles = theme => ({
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
    downloadButton:{
        textDecoration:'none',
        marginLeft:'20px',
        textTransform:'none',
    },
})

class QueryDendrogram extends React.Component {
    constructor(props) {
		super(props);
	}

	render(){
	    const { classes } = this.props;

	    if( this.props.dendrogram != undefined ){
	        return(
	            <div>
                    { this.props.dendrogram.png_file != undefined ?
                        <div>
                            <div style={{ margin: '20px', }} className={classes.displayCenter}>
                                <img src={this.props.dendrogram.svg_file} />
                            </div>
                            <br />
                            <div className={classes.displayCenter}>
                                <font>Download</font>
                            </div>
                            <br />
                            <div className={classes.displayCenter}>
                                <a download href={this.props.dendrogram.png_file} className={classes.downloadButton}>
                                    <Button variant="contained" color="default">Png</Button>
                                </a>
                                <a download href={this.props.dendrogram.pdf_file} className={classes.downloadButton}>
                                    <Button variant="contained" color="default">Pdf</Button>
                                </a>
                                <a download href={this.props.dendrogram.svg_file} className={classes.downloadButton}>
                                    <Button variant="contained" color="default">Svg</Button>
                                </a>
                                <a download href={this.props.dendrogram.newick_file} className={classes.downloadButton}>
                                    <Button variant="contained" color="default">Newick</Button>
                                </a>
                            </div>
                        </div> :
                        <div className={classes.displayCenter}>
                            <font style={{ margin: '50px 0 50px 0', fontSize: '20px', color: 'red'}}>
                                Data not found. Please input correct ID or try again later.
                            </font>
                        </div>
                    }
                </div>
	        )
	    }else{
	        return(
	            <div />
	        )
	    }
	}
}

export default withStyles(styles)(QueryDendrogram)