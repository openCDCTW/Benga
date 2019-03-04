import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
//Option's outlook
import PropTypes from 'prop-types';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
//
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import green from '@material-ui/core/colors/green';


const styles = theme => ({
    buttoncss:{
        color: theme.palette.getContrastText(green[600]),
        backgroundColor: green[500],
        '&:hover': {
            backgroundColor:green[600],
        },
    },
    selectcss: {
	    display: 'flex',
	    flexWrap: 'wrap',
	},
	formControl: {
	    margin: theme.spacing.unit,
	    minWidth: 240,
	},
	selectEmpty: {
	    marginTop: theme.spacing.unit * 2,
	},
    
})


class Tracking extends React.Component {

	constructor(props) {

		super(props);
		this.state = {
			to: "/tracking",
            upload_confirm: false,
            allele_db:"Vibrio_cholerae",
            profile_db:"Vibrio_cholerae",
		};

		this.djsConfig = {
            dictDefaultMessage:"Drop a cgMLST profile here",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 1,
            init:function(){
                this.on("addedfile", function(on_load_header_data){
                    // var fileSizeElement = on_load_header_data.previewElement.querySelector(".dz-size");
                    // fileSizeElement.parentNode.removeChild(fileSizeElement);
                });
                this.on("success", function(file,response){
                    file._removeLink.remove();
                    delete file._removeLink;
                    window.trackingID = response.id;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/tracking/sequence/'
        };

        this.dropzone = null;
	}

    handlePost() {
    
	    let fileCheck = this.dropzone.files.length;

	    if(fileCheck != 1){
	        alert('Please upload only 1 file');
	        return ;
	    }

	    this.dropzone.processQueue();
	    this.setState(state => ({ 
	    	upload_confirm: true,
	    	to: '/tracking_result' }));

	}

	submit(){

	    if (this.state.upload_confirm == false){
	        alert('Please upload a file first!');
	        return ;
	    }

	    var scheme = {};
	    scheme.id = window.trackingID ;
	    scheme.allele_db = this.state.allele_db;
        scheme.profile_db = this.state.profile_db;
	    fetch('api/tracking/tracking/', {
	        method:'POST',
	        headers: new Headers({'content-type': 'application/json'}),
	        body: JSON.stringify(scheme)
	    })
        .then(response => response.json());

        window.tabSwitch = true;
	}

	remove(){
	    this.dropzone.removeAllFiles();
	    this.setState(state => ({ to: '/tracking', upload_confirm: false }));

	}

	select_handleChange(event){
        if( event.target.value == 'Vibrio_cholerae'){
            this.setState(state => ({ 
                [event.target.name]: event.target.value,
                profile_db:"vibrio-profiles",
            }));
        }
    }

	render() {

		const config = this.componentConfig;
		const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        return (
            <div>
                <br />
                <br />
                <div>
                    <div style={{ float:'left', marginLeft:'10px' }}>
                        <form className={classes.selectcss} autoComplete="off">
                            <FormControl required className={classes.formControl} >
                              <InputLabel htmlFor="database-required">Database</InputLabel>
                                <Select
                                  value={this.state.allele_db}
                                  onChange={this.select_handleChange.bind(this)}
                                  name="allele_db"
                                  inputProps={{
                                    id: 'database-required',
                                  }}
                                  className={classes.selectEmpty}
                                  >
                                  <MenuItem value={'Vibrio_cholerae'}>Vibrio cholerae</MenuItem>
                                </Select>
                              <FormHelperText>Required</FormHelperText>
                            </FormControl>
                        </form>
                    </div>
                    <div style={{ float:'right', marginTop:'35px', marginRight:'25px' }}>
                        <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                                Remove all files
                                &nbsp;&nbsp;
                                <DeleteIcon />
                        </Button>
                    </div>
                </div>
                <br />
                <br />
                <br />
                <br />
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                </div>
                <br />
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" className ={classes.buttoncss} 
                     onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            Tracking
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <br />
                <br />
                <br />
            </div>
        );
    }
}

export default withStyles(styles)(Tracking);