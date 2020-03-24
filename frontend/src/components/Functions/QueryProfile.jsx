import React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import green from '@material-ui/core/colors/green';

const styles = theme => ({
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
    downloadButton: {
      color: theme.palette.getContrastText(green[900]),
      backgroundColor: green[600],
      '&:hover': {
        backgroundColor:green[900],
      },
      textTransform:'none',
    },

})

class QueryProfile extends React.Component {
    constructor(props) {
		super(props);
	}
	render(){
	    const { classes } = this.props;

        if( this.props.profile != undefined ){
            return (
                <div>
                    { this.props.profile.zip != undefined ?
                        <div className={classes.displayCenter}>
                            <a download href={this.props.profile.zip}
                             style={{ textDecoration:'none', margin: '70px' }}>
                                <Button variant="contained" className={classes.downloadButton} size="large">
                                    Download profiles (.zip)
                                    &nbsp;&nbsp;
                                    <DownloadIcon />
                                </Button>
                            </a>
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
            return (
                <div />
            )
        }

	}
}

export default withStyles(styles)(QueryProfile)