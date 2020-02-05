import React from 'react';
import ReactDOM from 'react-dom';
import Typography from '@material-ui/core/Typography';
import { withStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import green from '@material-ui/core/colors/green';
import VC from './DatabaseInfo_VC.jsx';
import NM from './DatabaseInfo_NM.jsx';

const styles = theme => ({
        downloadButton:{
          color: theme.palette.getContrastText(green[900]),
          backgroundColor: green[600],
          '&:hover': {
            backgroundColor:green[900],
          },
        },
});

class DatabaseInfo extends React.Component {

	constructor(props) {
		super(props);
	}

	render(){
		const { classes } = this.props;

		if(this.props.database=='VC'){
			return (<VC />)
		}else if(this.props.database=='NM'){
			return (<NM />)
		}
	}

}
export default withStyles(styles)(DatabaseInfo);
