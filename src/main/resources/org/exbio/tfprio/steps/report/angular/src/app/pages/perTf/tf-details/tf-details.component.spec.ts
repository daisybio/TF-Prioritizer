import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TfDetailsComponent } from './tf-details.component';

describe('TfDetailsComponent', () => {
  let component: TfDetailsComponent;
  let fixture: ComponentFixture<TfDetailsComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TfDetailsComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(TfDetailsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
