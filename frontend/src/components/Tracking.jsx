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
import blue from '@material-ui/core/colors/blue';
//search
import TextField from '@material-ui/core/TextField';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import SearchIcon from '@material-ui/icons/Search';
import Typography from '@material-ui/core/Typography';
import OutlinedInput from '@material-ui/core/OutlinedInput';

const styles = theme => ({
    buttoncss:{
        color: theme.palette.getContrastText(blue[600]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
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
        textField:{
        marginLeft: '19px',
        marginTop: '25px',
        width:'22.5%'
    },
    container:{
        display:'flex',
        flexWrap:'wrap',
    },
    numberField:{
        marginLeft: '19px',
        marginTop: '20px',
        marginRight: '10px',
        width:'110px'
    },
    countrySelect:{
        marginLeft: '19px',
        marginTop: '25px',
        width: '22.5%',
    }
    
})


class Tracking extends React.Component {

	constructor(props) {

		super(props);

        let nowYear = new Date();
        window.profile_db = "Vibrio_cholerae";

		this.state = {
            allele_db:"Vibrio_cholerae",
            yearError: false,
            country: '',
            labelWidth: 0,
            nowYear: nowYear.getFullYear(),
		};

		this.djsConfig = {
            dictDefaultMessage:"Drag a cgMLST profile here",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 1,
            maxFiles: 1,
            timeout: 0,
            init:function(){
                this.on("addedfile", function(file){
                    if(file.size <100000){
                        this.removeFile(file);
                    }
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("profile_db", window.profile_db);
                });
                this.on("success", function(file,response){
                    file._removeLink.remove();
                    delete file._removeLink;
                    window.trackingID = response.id;
                });
                this.on("maxfilesexceeded", function(file){
                    this.removeAllFiles();
                    this.addFile(file);
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/tracking/profile/'
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
        this.props.history.push("/cgMLST/tracking_result")
	}

	select_handleChange(event){
        if( event.target.value == 'Vibrio_cholerae'){
            this.setState(state => ({ 
                [event.target.name]: event.target.value
            }));
            window.profile_db = "Vibrio_cholerae";
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
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <a download href='https://drive.google.com/uc?export=download&id=1OgrY6Lt-93-gaYFEl7WbbUt6VLZ6eaC4' style={{ textDecoration:'none' }}>
                        <Button style={{ textTransform:'none' }} variant="contained" color="default">
                            Download &nbsp; example &nbsp; file
                        </Button>
                    </a>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Button variant="contained" className ={classes.buttoncss} 
                     onClick={this.handlePost.bind(this)}>
                        Submit
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                   <a download href='https://drive.google.com/uc?export=download&id=1_eQgWuGZxJAtOWxEgvdINyAJmyuENrza' style={{ textDecoration:'none' }}>
			<Button style={{ textTransform:'none' }} variant="contained" color="default">
                            Download &nbsp; cgMLST_tree_5048.pdf
			</Button>
                   </a>
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                   <a download href='https://drive.google.com/uc?export=download&id=1WwCoeIMBekda2z9BX1MMDqzW9y_pPfw7' style={{ textDecoration:'none' }}>
                        <Button style={{ textTransform:'none' }} variant="contained" color="default">
                            Download &nbsp; V.cholerae_genomes_list.zip
                        </Button>
                   </a>
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
